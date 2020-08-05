Hui-Syuan Yeh - 2571618
Akhmajon Makhsadov - 2571592
Shubhendu Jena - 2571723

Assignment H12

Problem 3 (Optic Flow Method of Lucas and Kanade)

Ans) We notice that for low values of integration scale {0.1,0.2} for both pig and sphere scenario, flow classification largely points to the aperture problem/normal flow scenario or to the scenario of having no flow information. For bigger values of integration scale [0.3,10], the flow classification improves and we are able to get flow field in almost all regions. However, with higher integration scale, we essentially convolve with gaussians of higher variances which leads to blurring of edges in flow fields.

Problem 4 (Variational Optic Flow Estimation)

Ans) We know that alpha which is the regularization parameter for the smoothness term controls the strength of the smoothness term in Horn-Schunk (HS) energy functional. We notice that as we increase value of alpha as alpha = {10,100,200,300,1000}, our penalty on the smoothness term increases and hence gradient of optic flow components (u,v) are forced to be small, hence making the flow smoother. We get the smoothest flow fields for alpha around 1000. Now, coming to the number of iterations, small number of iterations also don't give us good results as we need a decent amount of iterations for satisfactory propagation of flow to regions where optical flow computation is difficult due to low gradient value which renders the data term of the HS energy functional useless. We find that for alpha = 1000, we get not so good filling-in effect for number of iterations less than 400. We start getting better filling-in for number of iterations = 500 and above, with the most optimum filling-in occuring around number of iterations = {900,1000}
