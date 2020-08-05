Hui-Syuan Yeh - 2571618
Shubhendu Jena - 2571723
Akhmajon Makhsadov - 2571592


Assignment H7 (3)
Problem 3 (Edge and Corner detection)

a) 

Parameters for objects.pgm
sigma = 1
T1 = 1
T2 = 10

Parameters for objects are more flexible. Since the data is not noisy it is not necessary to even do any strong smoothing. Even without smooting we can safely compute derivatives. Another difference is that since the image is binary (black and white) it easier to find suitable T1 and T2 parameters. In other words we can just set T1 to 0 and T2 to any positive number to start tracing and we still get decent edges. And vice-versa, we can safely set T1 and T2 to pretty high values.

Parameters for pruebab1.pgm
sigma = 4
T1 = 4.5
T2 = 5

In this case we have to play with parameters more to get closed contours without artifacts.
Because of the strong noise we set sigma to 4. This value removes the noise and still preserves corners of the rectangle and the triangle.
T1 less than 4 for sigma 4 still leaks false edges, while T1 over 5 for sigma 4 kills true edges. We found that T1 = 4.5 works better than other values for choosen sigma. T2 we set to little higher value.


b) Corner detection

Parameters for tree.pgm
sigma = 1
rho = 0.25

Parameters for acros.pgm
sigma = 6
rho = 0.25

Parameters for stairs.pgm
sigma = 2
rho = 0.25


Problem 4: (Morphological Operations)

b)

To remove windows and doors in house.pgm:
Opening with m = 14

To create owl at night, owl.pgm:
Opening with m = 2 (the idea is that in the night you can only see owl's glowing eyes)

To separate image structures:
Angiogram.pgm - BTH with m=2
Fabric.pgm - WTH with m=2