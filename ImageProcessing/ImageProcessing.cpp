#include <iostream>
#include <Eigen/Dense>
#include <atlimage.h>
#include <algorithm>
#include <cmath>
#include <gurobi_c++.h>
#include <fftw3.h>

// X is Eigen::MatrixXd with rows = height, cols = width
Eigen::MatrixXd DctConvert(const Eigen::MatrixXd& input) {
    int rows = input.rows();
    int cols = input.cols();

    // FFTW real buffers (row‑major)
    double* in = fftw_alloc_real(rows * cols);
    double* out = fftw_alloc_real(rows * cols);

    // Copy Eigen data into FFTW buffer
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            in[i * cols + j] = input(i, j);             
            // Say, (1,2,3,4) is a matrix and 1,2 form the first row. 
            // The number of indices per row(which is also the number of columns) times the number of rows we want to cycle through to reach the x'th row(would be x - 1 if the index started at 1, but x since it started at 0)
            // lets us cycle to a new row in the first position, and adding to that the column index(if the first index started at 1 we would be adding y - 1) gives us what's stored in the (x, y) index

    // Create forward DCT‑II plan (2D)
    fftw_plan forward_plan = fftw_plan_r2r_2d(
        rows, cols,
        in, out,
        FFTW_REDFT10, FFTW_REDFT10,
        FFTW_ESTIMATE
    );

    // Execute forward DCT
    fftw_execute(forward_plan);

    // Copy result into Eigen
    Eigen::MatrixXd dct2_mat(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            dct2_mat(i, j) = out[i * cols + j];

    fftw_destroy_plan(forward_plan);
    fftw_free(in);
    fftw_free(out);

    return dct2_mat;
}

Eigen::MatrixXd DctInvert(const Eigen::MatrixXd& input) {
    int rows = input.rows();
    int cols = input.cols();

    // Allocate FFTW buffers
    double* in = fftw_alloc_real(rows * cols);
    double* out = fftw_alloc_real(rows * cols);

    // Copy input matrix into FFTW buffer (row-major)
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            out[i * cols + j] = input(i, j);

    // Create inverse DCT-III plan
    fftw_plan inverse_plan = fftw_plan_r2r_2d(
        rows, cols,
        out, in,
        FFTW_REDFT01, FFTW_REDFT01,
        FFTW_ESTIMATE
    );

    // Execute inverse DCT
    fftw_execute(inverse_plan);

    // Copy result into Eigen matrix
    Eigen::MatrixXd dct3_mat(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            dct3_mat(i, j) = in[i * cols + j];

    // Cleanup
    fftw_destroy_plan(inverse_plan);
    fftw_free(in);
    fftw_free(out);

    return dct3_mat;
}

Eigen::MatrixXd ZeroFunction(const Eigen::MatrixXd& input) {
    double lowbound = 1e-1;
    int rows = input.rows();
    int cols = input.cols();
    Eigen::MatrixXd output = input;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (input(i, j) < lowbound) { output(i, j) = 0.00; }
            }
        }
        return output;
} 
// zeroes out all values in indices lower than the bound

bool SaveMatrixAsImage(
    const Eigen::MatrixXd& X_R, 
    const Eigen::MatrixXd& X_G, 
    const Eigen::MatrixXd& X_B, 
    const wchar_t* outPath)
{
    int h = X_R.rows();
    int w = X_R.cols();

    CImage out; 
    out.Create(w, h, 24); // Creates buffer for the CImage object

    unsigned char* base = reinterpret_cast<unsigned char*>(out.GetBits()); // Gets a pointer to the buffer
    int pitch = out.GetPitch(); // bytes per row
    int bytesPerPixel = 3;

    for (int i = 0; i < h; ++i) {
        unsigned char* row = base + i * pitch;
        
        for (int j = 0; j < w; ++j) {
            double v_b = X_B(i, j);
            double clamped_v_b = std::clamp(v_b, 0.00, 255.00); // Not letting the values go outside of range
            unsigned char g_b = static_cast<unsigned char>(std::round(clamped_v_b)); // Rounds up the channel color value to store into the corresponding buffer 
            
            double v_g = X_G(i, j);
            double clamped_v_g = std::clamp(v_g, 0.00, 255.00);
            unsigned char g_g = static_cast<unsigned char>(std::round(clamped_v_g));
            
            double v_r = X_R(i, j);
            double clamped_v_r = std::clamp(v_r, 0.00, 255.00);
            unsigned char g_r = static_cast<unsigned char>(std::round(clamped_v_r));
            

            
            int offs = j * bytesPerPixel; 
            
            row[offs + 0] = g_b; // B
            row[offs + 1] = g_g; // G
            row[offs + 2] = g_r; // R
        }
    }

    HRESULT hr = out.Save(outPath); 
    return SUCCEEDED(hr);
}

Eigen::MatrixXd LPOptimize( 
    const Eigen::MatrixXd& B, 
    const Eigen::MatrixXd& C)
{
    int m = B.cols();
    int n = B.rows();
    int p = C.rows();

    GRBEnv env = GRBEnv(true);
    env.start();
    GRBModel model = GRBModel(env);

    std::vector<std::vector<GRBVar>> X(p, std::vector<GRBVar>(n));

    // Adds GRB variables according to the matrix
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < n; ++j) {
            X[i][j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }
    }

    model.update();


    for (int i = 0; i < p; ++i) {
        for (int k = 0; k < m; ++k) {
            GRBLinExpr sumExpr = 0;
            for (int j = 0; j < n; ++j) {
                sumExpr += X[i][j] * B(j, k);
            }
            model.addConstr(sumExpr == C(i, k));
        }
    }

    std::vector<std::vector<GRBVar>> absX(p, std::vector<GRBVar>(n));

    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < n; ++j) {
            absX[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

            model.addConstr(absX[i][j] >= X[i][j]);
            model.addConstr(absX[i][j] >= -X[i][j]);
        }
    }

    GRBLinExpr obj = 0;
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < n; ++j)
            obj += absX[i][j];

    model.setObjective(obj, GRB_MINIMIZE); 

    model.optimize();
    Eigen::MatrixXd A_recon(p, n);
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < n; ++j)
            A_recon(i, j) = X[i][j].get(GRB_DoubleAttr_X);
    return A_recon;
}

int main()
{
    CImage img;

    if (std::remove("Reconstruction/reconstructed.jpg") == 0) {
        std::cout << "Old image deleted.\n";
    }
    else {
        std::perror("Failed to delete file");
    }

    if (std::remove("Reconstruction/model.jpg") == 0) {
        std::cout << "Old model image deleted.\n";
    }
    else {
        std::perror("Failed to delete model file");
    }

    HRESULT hr = img.Load(L"Images/clay-banks-BnDI_MVomAI-unsplash.jpg");
    if (FAILED(hr)) {
        std::cerr << "Failed to load image!\n";
        return 1; // or handle error
    }
    else {
        std::cout << "Image loaded successfully.\n";
    }

    int m = 40;
    int n = img.GetWidth();
    int p = img.GetHeight();
    int bpp = img.GetBPP();           // bits per pixel
    int bytesPerPixel = bpp / 8;      // e.g., 1 for 8bpp, 3 for 24bpp, 4 for 32bpp
    int pitch = img.GetPitch();       // bytes per row (may be negative)
    unsigned char* base = reinterpret_cast<unsigned char*>(img.GetBits());

    // Define the true matrix
    Eigen::MatrixXd A_R(p, n);
    Eigen::MatrixXd A_G(p, n);
    Eigen::MatrixXd A_B(p, n);

    for (int i = 0; i < p; ++i) {
        unsigned char* row = (pitch < 0)
            ? base + i * pitch
            : base + (p - 1 - i) * (-pitch);
        for (int j = 0; j < n; ++j) {
            // pointer to the pixel's first byte
            unsigned char* px = row + j * bytesPerPixel;
            A_B(i, j) = static_cast<double>(px[0]);
            A_G(i, j) = static_cast<double>(px[1]);
            A_R(i, j) = static_cast<double>(px[2]);
        }
    }

    // Define a sampling matrix
    Eigen::MatrixXd B(n, m);

    // Generate random values 
    B = Eigen::MatrixXd::Random(n, m);

    Eigen::MatrixXd A_B_conv = ZeroFunction(DctConvert(A_B)); // Removes data that does not contribute much to the image
    Eigen::MatrixXd A_G_conv = ZeroFunction(DctConvert(A_G));
    Eigen::MatrixXd A_R_conv = ZeroFunction(DctConvert(A_R));

    Eigen::MatrixXd C_B = A_B_conv * B;
    Eigen::MatrixXd C_G = A_G_conv * B;
    Eigen::MatrixXd C_R = A_R_conv * B;

    Eigen::MatrixXd X_R_LP = LPOptimize(B, C_R);
    Eigen::MatrixXd X_G_LP = LPOptimize(B, C_G);
    Eigen::MatrixXd X_B_LP = LPOptimize(B, C_B);

    Eigen::MatrixXd A_B_inver = DctInvert(X_B_LP);
    Eigen::MatrixXd A_G_inver = DctInvert(X_G_LP);
    Eigen::MatrixXd A_R_inver = DctInvert(X_R_LP);

    A_B_inver /= 4 * p * n; // Normalizes values
    A_G_inver /= 4 * p * n;
    A_R_inver /= 4 * p * n;


    bool boolean2 = SaveMatrixAsImage(A_R_inver, A_G_inver, A_B_inver, L"Reconstruction/reconstructed.jpg");
    
    if (boolean2) { std::cout << "succeeded"; }
        else { std::cout << "failed"; }
 
    return 0;
}