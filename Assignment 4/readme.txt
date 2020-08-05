Hui-Syuan Yeh - 2571618
Shubhendu Jena - 2571723
Akhmajon Makhsadov - 2571592


Assignment H4 (4)
Problem 4 (Discrete Cosine Transform)

b)
1. On the DCT of 8x8 blocks image we can see rough structure of the original image. The source of this structure is high frequencies of the original image (higher coefs for higher frequencies fills 8x8 block with more gray values, while monotone blocks looks empty). The DCT of the whole image doesn't tell us much about its source.
2. Runtime for DCT of the whole image is way slower than DCT performed on the 8x8 blocks separately (30 sec vs 1 sec in our case). The reason is that the amount of basis vectors is increasing drastically along with computation cost for larger loops.

c)
1. We clearly see that high frequencies were removed from both spectra images.
2. On the resulting image produced by Menu 3 (DCT/IDCT of the whole image with removal of frequencies) we can see strong ringing effect produced by band pass filter. In spite that the hard edges are visually preserved better on this image, the overall visual quality is lower.
3. The image produced by Menu 4 (DCT/IDCT of 8x8 blocks with removal of frequencies) visually looks less damaged. But it is worth to mention, that if we zoom the image in, we see strong unpleasing borders between 8x8 blocks, which harms the image too. This kind of compression, with separation on 8x8 is not good enough to preserve hard, especially diagonal sharp edges across the image.

d)
1. Spectra of both images looks almost similar. But on the DCT with JPEG quantization we can see more low/middle frequencies preserved, and more high frequencies removed. In other words JPEG quantization is more focused to preserve lower frequencies.
2. Visual quality of the image produced by both images is better than images from previous options (Menu 3,4). In this case sharp edges are present.
3. IDCT with JPEG quantization looks better. More smart quantization strategy preserved those frequencies, which human's eye can see better. And sacrificed those frequencies to which human visual system is more tolerant.
