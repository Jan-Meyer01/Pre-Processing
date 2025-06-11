%-----------------------------------------------------------------------
% Job saved on 26-Feb-2025 15:15:21 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.hmri.denoise.subj.output.outdir = '<UNDEFINED>';
matlabbatch{1}.spm.tools.hmri.denoise.subj.pdw.mag_img = '<UNDEFINED>';
matlabbatch{1}.spm.tools.hmri.denoise.subj.pdw.phase_img = '';
matlabbatch{1}.spm.tools.hmri.denoise.subj.t1w.mag_img = '<UNDEFINED>';
matlabbatch{1}.spm.tools.hmri.denoise.subj.t1w.phase_img = '';
matlabbatch{1}.spm.tools.hmri.denoise.subj.mtw.mag_img = '<UNDEFINED>';
matlabbatch{1}.spm.tools.hmri.denoise.subj.mtw.phase_img = '';
matlabbatch{1}.spm.tools.hmri.denoise.subj.denoisingtype.lcpca_denoise.DNparameters.DNmetadata = 'yes';
matlabbatch{1}.spm.tools.hmri.denoise.subj.denoisingtype.lcpca_denoise.std = 1.05;
matlabbatch{1}.spm.tools.hmri.denoise.subj.denoisingtype.lcpca_denoise.ngbsize = 2;
matlabbatch{1}.spm.tools.hmri.denoise.subj.popup = false;
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = '<UNDEFINED>';
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.b1input = '<UNDEFINED>';
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.scafac = '<UNDEFINED>';
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.b1parameters.b1metadata = 'yes';
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT(1) = cfg_dep('Denoising: Denoised_mtw_magnitude', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','subj', '()',{1}, '.','DenoisedMagnitudemtw', '()',{':'}));
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD(1) = cfg_dep('Denoising: Denoised_pdw_magnitude', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','subj', '()',{1}, '.','DenoisedMagnitudepdw', '()',{':'}));
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1(1) = cfg_dep('Denoising: Denoised_t1w_magnitude', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','subj', '()',{1}, '.','DenoisedMagnitudet1w', '()',{':'}));
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;
