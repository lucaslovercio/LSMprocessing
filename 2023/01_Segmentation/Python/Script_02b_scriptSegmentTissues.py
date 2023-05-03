from TissueSegmentation.predict import segment_folder

# Model for 256 size tiles
model_dir       = '...h5'


frame_path      = ''
outputs_dir     = ''

segment_folder(model_dir, frame_path, outputs_dir)
