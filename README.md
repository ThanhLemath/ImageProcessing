Description
This is a project attempting to work with image optimization using GUROBI to perform minimization with the L1 norm, and DCT transformation using the three RGB channels
First, we convert a bitmap image into 3 matrices in three color channels
Then, we use apply forward Discrete Cosine Transform on them to induce sparsity in the matrix while still retaining meaningfull details
We multiply these matrices with another measurement matrix to create a system in this form
AB = C
Which is naturally underdetermined if B is not full rank
Then we apply the L1 norm to find the X solution matrices approximating A through GUROBI, leveraging the created sparsity
Then we perform backwards DCT, renormalize the image values and create the approximate output image

Dependencies
GUROBI (license), FFTW, ATL and Eigen are necessary to run the program
Install ATL separately
Put the library files for GUROBI, FFTW, and Eigen in the working directory and follow the .vcxproj file as the include and library files are stored deeper