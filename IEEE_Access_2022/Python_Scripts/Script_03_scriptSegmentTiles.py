from predict import segment_folder

model_dir = '' #H5 file path
frame_path = '' #Folder with tiles to segment tissues
outputs_dir = '' #Folder to save the segmented tiles

segment_folder(model_dir, frame_path, outputs_dir)
