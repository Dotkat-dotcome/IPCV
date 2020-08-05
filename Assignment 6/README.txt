Hui-Syuan Yeh - 2571618
Shubhendu Jena - 2571723
Akhmajon Makhsadov - 2571592


Assignment H6 (3)
Problem 3 (Linear Filters)

b)

LEOPARD
Noise is high frequency part of the image, this is why we have to get rid of it, by applying low pass filter.

To remove the noise from leopard image we used lowpass filter with sigma == 2.
With sigma == 1 some noise still remains, while sigma == 3 from out perspective oversmoothes the image.


ANGIOGRAM
To remove low-frequency background and focus on high-frequencies structures (vessels) we apply high-pass filter.

Sigma == 1 already gives good visible result. To get more prominent accent we applied Sigma == 2
Angiogram at sigma == 3 looks a bit oversharpened and has artifacts on boundaries which looks something like a ringing effect or halo.


TILE
According to the task we have to isolate the dark gaps.
It means we do not need any grain details of the tile (high-frequencies). Also we don't need any large scale variance (low-frequencies). In this case we are interested in certain frequency band.

We applied band filter with following parameters:
Sigma1 == 2.5 (enough to keep the tile structure)
Sigma2 == 2.0 (enough to remove grain)