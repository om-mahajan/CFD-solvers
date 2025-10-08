// ME22B174. Install EIGEN and run. I have given a compiled version, run using ./out1.exe
#include "SIMPLE.h"

double SIMPLE::calculateStep() {
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N - 1; j++) {

            if (checkBoundaries(i, j) == 1) {
                continue;
            }

            double Fe = 0.5 * rho * (u(i, j) + u(i, j + 1)) * hx;
            double Fw = 0.5 * rho * (u(i, j) + u(i, j - 1)) * hx;
            double Fn = 0.5 * rho * (v(i - 1, j) + v(i - 1, j + 1)) * hy;
            double Fs = 0.5 * rho * (v(i, j) + v(i, j + 1)) * hy;

            double Pe = Fe / D;

            if (abs(Pe) > 1) {
                //std::cout << Pe << std::endl;
            }
            //double aE = -0.5 * Fe + D;
            //double aW = 0.5 * Fw + D;
            //double aN = -0.5 * Fn + D;
            //double aS = 0.5 * Fs + D;

            //double aE = D + takeMax(0, -Fe);
            //double aW = D + takeMax(0, Fw); 
            //double aN = D + takeMax(0, -Fn); 
            //double aS = D + takeMax(0, Fs); 

            double aE = takeMax(-Fe, D - 0.5 * Fe, 0);
            double aW = takeMax(Fw, D + 0.5 * Fw, 0);
            double aN = takeMax(-Fn, D - 0.5 * Fn, 0);
            double aS = takeMax(Fs, D + 0.5 * Fs, 0);

            double aP = (aW + aE + aN + aS + Fe - Fw + Fn - Fs);

            double Ae = -hx;

            dE(i, j) = Ae / aP;

            uStar(i, j) = (aE * u(i, j + 1) + aW * u(i, j - 1) + aN * u(i - 1, j) + aS * u(i + 1, j)) / aP + dE(i, j) * (p(i, j + 1) - p(i, j));

        }
    }

    for (int i = 1; i < M - 1; i++) {
        for (int j = 1; j < N; j++) {

            if (checkBoundaries(i, j) == 1) {
                continue;
            }

            double Fe = 0.5 * rho * (u(i, j) + u(i + 1, j)) * hx;
            double Fw = 0.5 * rho * (u(i, j - 1) + u(i + 1, j - 1)) * hx;
            double Fn = 0.5 * rho * (v(i - 1, j) + v(i, j)) * hy;
            double Fs = 0.5 * rho * (v(i, j) + v(i + 1, j)) * hy;

            //double aE = -0.5 * Fe + D;
            //double aW = 0.5 * Fw + D;
            //double aN = -0.5 * Fn + D;
            //double aS = 0.5 * Fs + D;

            ////double aE = D + takeMax(0, -Fe);
            ////double aW = D + takeMax(0, Fw);
            ////double aN = D + takeMax(0, -Fn);
            ////double aS = D + takeMax(0, Fs);

            double aE = takeMax(-Fe, D - 0.5 * Fe, 0);
            double aW = takeMax(Fw, D + 0.5 * Fw, 0);
            double aN = takeMax(-Fn, D - 0.5 * Fn, 0);
            double aS = takeMax(Fs, D + 0.5 * Fs, 0);

            double aP = (aW + aE + aN + aS + Fe - Fw + Fn - Fs);

            double An = -hy;

            dN(i, j) = An / aP;

            vStar(i, j) = (aE * v(i, j + 1) + aW * v(i, j - 1) + aN * v(i - 1, j) + aS * v(i + 1, j)) / aP + dN(i, j) * (p(i, j) - p(i + 1, j));
        }
    }

    setVelocityBoundaryConditions(uStar, vStar);

    pStar = Eigen::MatrixXd::Zero(M + 1, N + 1);

    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {

            if (checkBoundaries(i, j) == 1) {
                continue;
            }

            double aE = -rho * dE(i, j) * hx;
            double aW = -rho * dE(i, j - 1) * hx;
            double aN = -rho * dN(i - 1, j) * hy;
            double aS = -rho * dN(i, j) * hy;

            double aP = aE + aW + aN + aS;

            b(i, j) = -rho * (uStar(i, j) - uStar(i, j - 1)) * hx + rho * (vStar(i, j) - vStar(i - 1, j)) * hy;


            pStar(i, j) = (aE * pStar(i, j + 1) + aW * pStar(i, j - 1) + aN * pStar(i - 1, j) + aS * pStar(i + 1, j) + b(i, j)) / aP;
        }
    }

    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {

            if (checkBoundaries(i, j) == 1) {
                continue;
            }

            p(i, j) = p(i, j) + pAlpha * pStar(i, j);
        }
    }

    setPressureBoundaryConditions(p);

    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N - 1; j++) {

            if (checkBoundaries(i, j) == 1) {
                continue;
            }

            u(i, j) = uStar(i, j) + uvAlpha * dE(i, j) * (pStar(i, j + 1) - pStar(i, j));
        }
    }

    for (int i = 1; i < M - 1; i++) {
        for (int j = 1; j < N; j++) {

            if (checkBoundaries(i, j) == 1) {
                continue;
            }

            v(i, j) = vStar(i, j) + uvAlpha * dN(i, j) * (pStar(i, j) - pStar(i + 1, j));
        }
    }

    setVelocityBoundaryConditions(u, v);

    double error = 0;

    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {

            if (checkBoundaries(i, j) == 1) {
                continue;
            }

            error += abs(b(i, j));
        }
    }

    return error;
}

double SIMPLE::takeMax(double A, double B) {
    double max = A; 
    if (B > max) {
        max = B;
    }
    return max;
}

double SIMPLE::takeMax(double A, double B, double C) {
    double max = A; 
    if (B > max) {
        max = B;
    }

    if (C > max) {
        max = C;
    }

    return max;
}

void SIMPLE::setVelocityBoundaryConditions(Eigen::MatrixXd& uIn, Eigen::MatrixXd& vIn) {

    for (int j = 0; j < N; j++) {
        uIn(0, j) = uIn(1, j); 
        vIn(0, j) = 0;
        vIn(M - 1, j) = 0;
    }

    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {

            if (checkBoundaries(i, j) == 1) {
                uIn(i, j) = 0;
                vIn(i, j) = 0;
            }
        }
    }

    for (int j = 0; j < N; j++) {
        uIn(0, j) = uIn(1, j);
        uIn(M, j) = uIn(M-1, j);
    }

    for (int i = 0; i < M; i++) {
        uIn(i, 0) = 10;
        uIn(i, N - 1) = uIn(i, N - 2);

        vIn(i, 0) = 0;
        vIn(i, N) = 0;
    }
}

void SIMPLE::setPressureBoundaryConditions(Eigen::MatrixXd& pIn) {

    for (int j = 1; j < N + 1; j++) {
        pIn(0, j) = pIn(1, j);
        pIn(M, j) = pIn(M - 1, j);
    }

    for (int i = 1; i < M + 1; i++) {
        pIn(i, 0) = pIn(i, 1);
        pIn(i, N) = pIn(i, N - 1);
    }
}

double SIMPLE::checkBoundaries(int i, int j) {
    // Cylinder centered at (cx, cy) with radius r
    double cx = M / 2.0;  // center i-index
    double cy = N / 2.0;  // center j-index
    double r = std::min(M, N) / 6.0;  // radius in grid units (adjust as needed)

    // Compute squared distance from (i, j) to (cx, cy)
    double dist2 = (i - cx) * (i - cx) + (j - cy) * (j - cy);

    // Return 1 if inside cylinder (obstacle), else 0
    return (dist2 <= r * r) ? 1.0 : 0.0;
}



void SIMPLE::saveMatrix(Eigen::MatrixXd inputMatrix, std::string fileName) {

    std::ofstream file;

    file.open(fileName + ".txt");

    //for (int i = 0; i < M; i++) {
    //    for (int j = 0; j < N; j++) {

    //        u(i, j) = 0.5 * (u(i, j) + u(i + 1, j));
    //        v(i, j) = 0.5 * (v(i, j) + v(i, j + 1));
    //    }
    //}

    //for (int i = 0; i < M; i++) {
    //    for (int j = 0; j < N; j++) {
    //        uvAbs(i, j) = sqrt(u(i, j) * u(i, j) + v(i, j) * v(i, j));
    //    }
    //}

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            file << "\t" << inputMatrix(i, j);
        }
        file << std::endl;
    }
}

void SIMPLE::saveAll() {

    Eigen::MatrixXd uOut = u; 
    Eigen::MatrixXd vOut = v;
    Eigen::MatrixXd pOut = p;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            uOut(i, j) = 0.5 * (u(i, j) + u(i + 1, j));
            vOut(i, j) = 0.5 * (v(i, j) + v(i, j + 1));
            pOut(i, j) = 0.25 * (p(i, j) + p(i, j + 1) + p(i + 1, j) + p(i + 1, j + 1));
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            pOut(i, j) = 0.25 * (p(i, j) + p(i, j + 1) + p(i + 1, j) + p(i + 1, j + 1));
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            uvAbs(i, j) = sqrt(u(i, j) * u(i, j) + v(i, j) * v(i, j));
        }
    }

    saveMatrix(uOut, "u");
    saveMatrix(vOut, "v");
    saveMatrix(pOut, "p");
    saveMatrix(uvAbs, "speed");

}

void SIMPLE::paintBoundaries() {

    Eigen::MatrixXd BCs = Eigen::MatrixXd::Zero(M, N);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            BCs(i, j) = checkBoundaries(i, j);
        }
    }

    saveMatrix(BCs, "BC");
}

SIMPLE::SIMPLE() {

    uStar = Eigen::MatrixXd::Zero(M + 1, N);
    vStar = Eigen::MatrixXd::Zero(M, N + 1);
    pStar = Eigen::MatrixXd::Zero(M + 1, N + 1);

    u = Eigen::MatrixXd::Zero(M + 1, N);
    v = Eigen::MatrixXd::Zero(M, N + 1);
    p = Eigen::MatrixXd::Zero(M + 1, N + 1);

    b = Eigen::MatrixXd::Zero(M + 1, N + 1);
    dE = Eigen::MatrixXd::Zero(M + 1, N);
    dN = Eigen::MatrixXd::Zero(M, N + 1);

    uvAbs = Eigen::MatrixXd::Zero(M + 1, N + 1);
    
    u(1, 1) = 1; 
    p(M, N) = 1;

}

void SIMPLE::runIterations() {
    int numberIterations = 0;
    double error = 1;
    paintBoundaries();
    std::ofstream diagnosticsFile; 
    diagnosticsFile.open("Convergence.txt"); 

    while (error > epsilon && numberIterations < maxIterations) {
        error = calculateStep();

        if (numberIterations % 1000 == 0) {
            std::cout << "Iteration " << numberIterations << "\t Error = " << error << std::endl;
            diagnosticsFile << numberIterations << "\t" << error << std::endl; 
        }

        if (numberIterations % 5000 == 0) {
            saveAll(); 
            
        }

        numberIterations += 1;
    }
    saveAll(); 
}


int main()
{

    SIMPLE CFD;
    CFD.runIterations();

}
