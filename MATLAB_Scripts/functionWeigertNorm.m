function slice_original_quantile = functionWeigertNorm(slice_original_quantile)
vQuantiles = quantile(slice_original_quantile(:),1000);
quantile999low = vQuantiles(1);%double(prctile(slice_original_quantile(:),0.001));
quantile999max = vQuantiles(1000);%double(prctile(slice_original_quantile(:),0.999));
slice_original_quantile(slice_original_quantile<quantile999low) = quantile999low;
slice_original_quantile(slice_original_quantile>quantile999max) = quantile999max;

if quantile999max > 10%to avoid empty image becoming white or enhacing noise
    
    slice_original_quantile = (double(slice_original_quantile - quantile999low) )./ double(quantile999max-quantile999low);
else
    slice_original_quantile = zeros(size(slice_original_quantile));
    
end

