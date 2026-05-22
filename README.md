Description<br>
This is a project attempting to work with image optimization using GUROBI to perform minimization with the L1 norm, and DCT transformation using the three RGB channels<br>
First, we convert a bitmap image into 3 matrices in three color channels<br>
Then, we use apply forward Discrete Cosine Transform on them to induce sparsity in the matrix while still retaining meaningfull details<br>
We multiply these matrices with another measurement matrix to create a system in the form<br>
AB = C
Which is naturally underdetermined if B is not full rank<br>
Then we apply the L1 norm to find the X solution matrices approximating A through GUROBI, leveraging the created sparsity<br>
Then we perform backwards DCT, renormalize the image values and create the approximate output image<br>

Dependencies<br>
GUROBI (license), FFTW, ATL and Eigen are necessary to run the program<br>
Install ATL separately<br>
Put the library files for GUROBI, FFTW, and Eigen in the working directory and follow the .vcxproj file as the include and library files are stored deeper<br>