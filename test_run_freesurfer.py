import shlex
import subprocess
import sys
from os.path import isfile, split, join
from os import access, X_OK, environ, getenv, pathsep

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

def run_cmd(cmd):
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
    progname=which(clist[0]) #join('~','freesurfer','bin','mri_synthstrip.sh')
    if (progname) is None:
        print('ERROR: '+ clist[0] +' not found in path!')
        sys.exit(1)
    clist[0]=progname
    cmd = ' '.join(clist)
    print('#@# Command: ' + cmd+'\n')

    args = shlex.split(cmd)
    process = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()  # Wait for the process to complete
    stdout, stderr = process.communicate()
    print(str(stdout))
    if stderr != b'':
        print('ERROR: '+ str(stderr))

    return stdout


save_name = '/home/janmeyer/Pre-Processing/processed_data/Test/sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii'
mask_name = '/home/janmeyer/Pre-Processing/processed_data/Test/brain_mask_sub-tle001.nii'

cmd = 'mri_synthstrip -i {} -o {} -m {}'.format(save_name,save_name,mask_name)
fs_output = run_cmd(cmd)