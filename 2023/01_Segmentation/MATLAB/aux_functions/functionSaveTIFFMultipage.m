function functionSaveTIFFMultipage(volume,filenameOutput,bitdepth)
[h,w,d]=size(volume);

if bitdepth == 8
    volume = uint8(volume);
else
    volume = uint16(volume);
end

imwrite(volume(:,:,1), filenameOutput);
for i=2:d
    imwrite(volume(:,:,i), filenameOutput, 'writemode', 'append');
end


end

