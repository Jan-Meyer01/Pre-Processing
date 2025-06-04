from os import walk, scandir, access, X_OK, environ, getenv, pathsep, mkdir
from os.path import join, isdir, islink, isfile, split
import numpy as np
from numpy import inf
import torch
import torch.nn as nn
import surfa as sf
import time
import shlex
import subprocess
import sys

## helper scripts for finding the nifit images in a given folder

def is_image_file(filename):
    '''
    Check whether image is actually a nifti file
    '''
    IMG_EXTENSIONS = ['.nii', '.nii.gz']        # only nifti files
    return any(filename.endswith(extension) for extension in IMG_EXTENSIONS)

def get_image_paths(dir, contrasts):
    '''
    Returns an array with the paths of all nifti images from a given directory
    '''
    images_paths = []
    assert isdir(dir) or islink(dir), '%s is not a valid directory' % dir

    # if subfolder shall be ignored
    fnames = [f.path for f in scandir(dir)]

    # if subfolders are wanted
    #for root, _, fnames in sorted(walk(dir, followlinks=True)):
    for fname in fnames:
        if is_image_file(fname):
            #path = join(root, fname)
            path = fname
            if type(contrasts) == type(None):
                images_paths.append(path)
            elif contrasts == []:  
                images_paths.append(path)  
            else:
                for contrast in contrasts:
                    if path.find(contrast) != -1:
                        images_paths.append(path)
    return images_paths

## helper functions for processing and converting travelling head R and T maps
def process_and_convert_R_to_T_travelling_head(R_image, brain_mask, cutoff):
    '''
    Function for processing and converting R maps to T maps for the travelling head dataset
    Expects an R map [np.array] and brain mask [np.array] with the same size as well as a cutoff for the thresholding
    Returns both the processed R map [np.array] and the converted T map [np.array] (same size as the inputs)
    '''
    # set negative and inf values to NaN
    R_image[R_image<0]    = np.NaN
    R_image[R_image==inf] = np.NaN
    # convert from R to T using a copy of the volume
    R_image_convert = R_image
    # thresholding for feasable values (we know which T values to expect in human tissue)
    R_image_convert[R_image_convert<cutoff] = np.NaN
    # invert the map
    R_image_convert = 1/R_image_convert
    # scale because values in ms are needed for JEMRIS simulation
    #volume = volume * 1000
    # remove NaN values after conversion
    R_image_convert[np.isnan(R_image_convert)] = 0
    # extract the brain from the images
    R_image[brain_mask==0] = 0
    R_image_convert[brain_mask==0] = 0

    return R_image, R_image_convert

def process_and_convert_T_to_R_travelling_head(T_image, brain_mask, cutoff):
    '''
    Function for processing and converting T maps to R maps for the travelling head dataset
    Expects an T map [np.array] and brain mask [np.array] with the same size as well as a cutoff for the thresholding
    Returns both the processed T map [np.array] and the converted R map [np.array] (same size as the inputs)
    '''
    # set negative and inf values to NaN
    T_image[T_image<0]    = np.NaN
    T_image[T_image==inf] = np.NaN
    # convert from T to R using a copy of the volume
    T_image_convert = T_image
    # thresholding for feasable values (we know which T values to expect in human tissue)
    T_image_convert[T_image_convert<cutoff] = np.NaN
    # invert the map
    T_image_convert = 1/T_image_convert
    # scale because values in ms are needed for JEMRIS simulation
    #volume = volume * 1000
    # remove NaN values after conversion
    T_image_convert[np.isnan(T_image_convert)] = 0
    # extract the brain from the images
    T_image[brain_mask==0] = 0
    T_image_convert[brain_mask==0] = 0

    return T_image, T_image_convert


## helper scripts for performing the image pre-processing

def process_volume_data(vol, debug=False):
    '''
    Given a 3D volume [np.array], remove all negative, NaN and inf values, 
    extract the brain and remove any remaining outliers in the data.
    If debug=True, all steps of the pre-processing are returned for debugging and visualization.
    '''
    
    print('Cleaning data')
    # clean data (remove negative, NaN and inf values)
    if debug:
        start = time.time()
        vol_clean = data_cleaning(vol)
        print('Cleaning took {} seconds'.format(time.time()-start))
    else:
        vol = data_cleaning(vol)
    
    print('Removing Outliers')
    # remove outliers
    if debug:
        start = time.time()
        vol_no_outliers = clip_outliers(vol_clean, threshold=3)
        print('Outlier removal took {} seconds'.format(time.time()-start))
    else:
        vol = clip_outliers(vol, threshold=3)        
    
    print('Skull-stripping')
    # skull stripping to extract the brain
    if debug:
        start = time.time()
        vol_skull_stripped = skull_stripping(vol_no_outliers, debug=debug)
        print('Skull-stripping took {} seconds'.format(time.time()-start))
    else:
        vol = skull_stripping(vol, debug=debug)
    
    if debug:
        return vol_clean, vol_skull_stripped, vol_no_outliers
    else: 
        return vol

def data_cleaning(vol, replacement_value=0):
    '''
    Set negative, NaN and inf values from a given 3D volume [np.array] to a replacement value [int or float].
    '''

    # zero out negative values (not quite 0, but 0.0001 in case of inversion)
    vol[vol < 0] = 0.0001 #replacement_value
    # NaN values
    vol[np.isnan(vol)] = replacement_value
    # and inf values
    vol[vol == inf] = replacement_value
    
    return vol

def clip_outliers(vol, threshold=0.9975):
    '''
    Clip outliers from a given 3D volume [np.array] based on the histogram and a given threshold [float] (in percent) of values after which outliers are clipped to the highest value. 
    '''

    # get histogram of the volume with 100.000 bins for better accuracy
    hist, bin_edges = np.histogram(np.ndarray.flatten(vol), bins=100000, density=False)
    index = 0
    cumultative_value = hist[index]
    
    # stop accumulating the bins once the threshold is reached
    while cumultative_value < np.sum(hist)*threshold:
        index += 1
        cumultative_value += hist[index]

    # and remove all values above the indexed bin
    vol[vol > bin_edges[index]] = bin_edges[index]

    return vol

def remove_outliers_median(vol, threshold=3):
    '''
    Remove outliers from a given 3D volume [np.array] based on the data distribution and a given threshold [int] for 
    the number of standard deviation away from the median value that are still excepted as normal values.
    '''

    # flatten the array to a single dimension, but remember the shape of the original volume
    vol_size = vol.shape
    vol = np.ndarray.flatten(vol)
    
    # set all values that fall outside the threshold times std to zero
    vol[abs(vol - np.median(vol)) > threshold * np.std(vol)] = 0

    # resize the processed volume back to the original size
    vol = np.reshape(vol,vol_size)

    return vol


## helper functions to call freesurfers SynthStrip model through the command line
def which(program):
    def is_exe(fpath):
        return isfile(fpath) and access(fpath, X_OK)

    fpath, _ = split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in environ["PATH"].split(pathsep):
            path = path.strip('"')
            exe_file = join(path, program)
            if is_exe(exe_file):
                return exe_file
        home = getenv('FREESURFER_HOME') #'/home/janmeyer/freesurfer'
        program_dir = join(home,'bin',program)
        if is_exe(program_dir):
            return program_dir
        program_dir = join('.',program)
        if is_exe(program_dir):
            return program_dir

    return None

def run_freesurfer(cmd):
    """
    execute the comand
    """
    environ["FREESURFER_HOME"] = "/home/janmeyer/freesurfer"
    command = ['/bin/sh', '/home/janmeyer/freesurfer/SetUpFreeSurfer.sh']
    
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    command = ['/bin/sh', '/home/janmeyer/freesurfer/FreeSurferEnv.sh']
    
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    clist = cmd.split()
    progname=which(clist[0]) 
    if (progname) is None:
        print('ERROR: '+ clist[0] +' not found in path!')
        sys.exit(1)
    clist[0]=progname
    cmd = ' '.join(clist)
    print('#@# Command: ' + cmd+'\n')

    args = shlex.split(cmd)
    process = subprocess.call(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #process.wait()  # Wait for the process to complete
    #stdout, stderr = process.communicate()
    print(str(stdout))
    if stderr != b'':
        print('ERROR: '+ str(stderr))

    return stdout

## helper function to use SynthStrip directly --> results look worse than the command line version?!?


def skull_stripping(vol, debug=False):
    '''
    SynthStrip for skull-stripping - Code taken from:
    ----------------------------------------------------
    SynthStrip: Skull-Stripping for Any Brain Image
    A Hoopes, JS Mora, AV Dalca, B Fischl, M Hoffmann
    NeuroImage 206 (2022), 119474
    https://doi.org/10.1016/j.neuroimage.2022.119474

    Website: https://synthstrip.io
    '''

    # necessary for speed gains (I think)
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.deterministic = True

    # configure device
    if torch.cuda.is_available():
        device = torch.device('cuda')
        if debug:
            print('running on GPU')
    else:
        device = torch.device('cpu')
        if debug:
            print('running on CPU')

    def extend_sdt(sdt, border=1):
        """Extend SynthStrip's narrow-band signed distance transform (SDT).

        Recompute the positive outer part of the SDT estimated by SynthStrip, for
        borders that likely exceed the 4-5 mm band. Keeps the negative inner part
        intact and only computes the outer part where needed to save time.

        Parameters
        ----------
        sdt : sf.Volume
            Narrow-band signed distance transform estimated by SynthStrip.
        border : float, optional
            Mask border threshold in millimeters.

        Returns
        -------
        sdt : sf.Volume
            Extended SDT.

        """
        if border < int(sdt.max()):
            return sdt

        # Find bounding box.
        mask = sdt < 1
        keep = np.nonzero(mask)
        low = np.min(keep, axis=-1)
        upp = np.max(keep, axis=-1)

        # Add requested border.
        gap = int(border + 0.5)
        low = (max(i - gap, 0) for i in low)
        upp = (min(i + gap, d - 1) for i, d in zip(upp, mask.shape))

        # Compute EDT within bounding box. Keep interior values.
        ind = tuple(slice(a, b + 1) for a, b in zip(low, upp))
        out = np.full_like(sdt, fill_value=100)
        out[ind] = sf.Volume(mask[ind]).distance()
        out[keep] = sdt[keep]

        return sdt.new(out)

    # configure model
    class StripModel(nn.Module):

        def __init__(self,
                    nb_features=16,
                    nb_levels=7,
                    feat_mult=2,
                    max_features=64,
                    nb_conv_per_level=2,
                    max_pool=2,
                    return_mask=False):

            super().__init__()

            # dimensionality
            ndims = 3

            # build feature list automatically
            if isinstance(nb_features, int):
                if nb_levels is None:
                    raise ValueError('must provide unet nb_levels if nb_features is an integer')
                feats = np.round(nb_features * feat_mult ** np.arange(nb_levels)).astype(int)
                feats = np.clip(feats, 1, max_features)
                nb_features = [
                    np.repeat(feats[:-1], nb_conv_per_level),
                    np.repeat(np.flip(feats), nb_conv_per_level)
                ]
            elif nb_levels is not None:
                raise ValueError('cannot use nb_levels if nb_features is not an integer')

            # extract any surplus (full resolution) decoder convolutions
            enc_nf, dec_nf = nb_features
            nb_dec_convs = len(enc_nf)
            final_convs = dec_nf[nb_dec_convs:]
            dec_nf = dec_nf[:nb_dec_convs]
            self.nb_levels = int(nb_dec_convs / nb_conv_per_level) + 1

            if isinstance(max_pool, int):
                max_pool = [max_pool] * self.nb_levels

            # cache downsampling / upsampling operations
            MaxPooling = getattr(nn, 'MaxPool%dd' % ndims)
            self.pooling = [MaxPooling(s) for s in max_pool]
            self.upsampling = [nn.Upsample(scale_factor=s, mode='nearest') for s in max_pool]

            # configure encoder (down-sampling path)
            prev_nf = 1
            encoder_nfs = [prev_nf]
            self.encoder = nn.ModuleList()
            for level in range(self.nb_levels - 1):
                convs = nn.ModuleList()
                for conv in range(nb_conv_per_level):
                    nf = enc_nf[level * nb_conv_per_level + conv]
                    convs.append(ConvBlock(ndims, prev_nf, nf))
                    prev_nf = nf
                self.encoder.append(convs)
                encoder_nfs.append(prev_nf)

            # configure decoder (up-sampling path)
            encoder_nfs = np.flip(encoder_nfs)
            self.decoder = nn.ModuleList()
            for level in range(self.nb_levels - 1):
                convs = nn.ModuleList()
                for conv in range(nb_conv_per_level):
                    nf = dec_nf[level * nb_conv_per_level + conv]
                    convs.append(ConvBlock(ndims, prev_nf, nf))
                    prev_nf = nf
                self.decoder.append(convs)
                if level < (self.nb_levels - 1):
                    prev_nf += encoder_nfs[level]

            # now we take care of any remaining convolutions
            self.remaining = nn.ModuleList()
            for num, nf in enumerate(final_convs):
                self.remaining.append(ConvBlock(ndims, prev_nf, nf))
                prev_nf = nf

            # final convolutions
            if return_mask:
                self.remaining.append(ConvBlock(ndims, prev_nf, 2, activation=None))
                self.remaining.append(nn.Softmax(dim=1))
            else:
                self.remaining.append(ConvBlock(ndims, prev_nf, 1, activation=None))

        def forward(self, x):

            # encoder forward pass
            x_history = [x]
            for level, convs in enumerate(self.encoder):
                for conv in convs:
                    x = conv(x)
                x_history.append(x)
                x = self.pooling[level](x)

            # decoder forward pass with upsampling and concatenation
            for level, convs in enumerate(self.decoder):
                for conv in convs:
                    x = conv(x)
                if level < (self.nb_levels - 1):
                    x = self.upsampling[level](x)
                    x = torch.cat([x, x_history.pop()], dim=1)

            # remaining convs at full resolution
            for conv in self.remaining:
                x = conv(x)

            return x

    class ConvBlock(nn.Module):
        """
        Specific convolutional block followed by leakyrelu for unet.
        """

        def __init__(self, ndims, in_channels, out_channels, stride=1, activation='leaky'):
            super().__init__()

            Conv = getattr(nn, 'Conv%dd' % ndims)
            self.conv = Conv(in_channels, out_channels, 3, stride, 1)
            if activation == 'leaky':
                self.activation = nn.LeakyReLU(0.2)
            elif activation == None:
                self.activation = None
            else:
                raise ValueError(f'Unknown activation: {activation}')

        def forward(self, x):
            out = self.conv(x)
            if self.activation is not None:
                out = self.activation(out)
            return out

    with torch.no_grad():
        model = StripModel()
        model.to(device)
        model.eval()

    # load model weights
    modelfile = join('/home','janmeyer','freesurfer', 'models', 'synthstrip.1.pt')
    #modelfile = join('.', 'models', 'synthstrip.1.pt')
    checkpoint = torch.load(modelfile, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])

    # convert np.array to surfa volume
    vol_sf = sf.Volume(vol)
    if debug:
        print(f'Processing frame (of {vol_sf.nframes}):', end=' ', flush=True)

    # loop over frames (try not to keep too much data in memory)
    dist = []
    mask = []
    for f in range(vol_sf.nframes):
        if debug:
            print(f + 1, end=' ', flush=True)
        frame = vol_sf.new(vol_sf.framed_data[..., f])

        # conform, fit to shape with factors of 64
        conformed = frame.conform(voxsize=1.0, dtype='float32', method='nearest', orientation='LIA')
        conformed = conformed.crop_to_bbox()
        target_shape = np.clip(np.ceil(np.array(conformed.shape[:3]) / 64).astype(int) * 64, 192, 320)
        conformed = conformed.reshape(target_shape)

        # normalize
        inp = torch.from_numpy(conformed.data).to(device).unsqueeze(0).unsqueeze(0)
        inp -= inp.min()
        inp /= inp.quantile(0.99)
        inp = inp.clamp(0, 1)

        # predict the sdt
        with torch.no_grad():
            sdt = model(inp).squeeze().cpu()

        # extend the sdt if needed, unconform
        border = 1
        sdt = extend_sdt(conformed.new(sdt), border=border)
        sdt = sdt.resample_like(vol_sf, fill=100)
        dist.append(sdt)

        # extract mask, find largest CC to be safe
        mask.append((sdt < border).connected_component_mask(k=1, fill=True))

    # combine frames and end line
    dist = sf.stack(dist)
    mask = sf.stack(mask)

    # mask the volume
    vol_tmp = vol
    vol_tmp[mask == 0] = np.min([0, vol_tmp.min()])
    if debug:
        print('Done!')

    return vol_tmp
