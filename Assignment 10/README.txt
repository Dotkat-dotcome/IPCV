Hui-Syuan Yeh - 2571618
Shubhendu Jena - 2571723
Akhmajon Makhsadov - 2571592

Assignment H10
Problem 3 (Texture Inpainting)
a) 
Parameter "m" is basically the patch size of the patch that we use to extract the center pixel of the most similar patch in the known region so as to inpaint the unknown region. If we keep the value of m high, we reduce the number of patch comparisons and our range of choices is limited. This is beneficial if our image is regular and the hence the unknown region is captured by comparing patches of bigger sizes. This is demonstrated by the case of circles.pgm as in part a). However, if as our image goes on becoming more and more irregular, we ought to decrease our patch size as the inpainting process would require us to capture local details better which would work well if our similarity patch sizes were smaller.   


Problem 4 (Deconvolution with Wiener Filtering)
b) K = 0.001
c) K = 0.00001
d) K = 0.005

The image size is not power of 2 ==> the algorithm cannot use FFT, and uses DFT instead. 
This is main cause of slowering the debluring.
