#include "fdtd-pml-2d.h"

void FdtdPml2D::SetParams(std::vector<std::vector<double>> &eps,
                          std::vector<std::vector<double>> &mu,
                          std::vector<std::vector<double>> &sigma,
                          size_t new_src_position_row,
                          size_t new_src_position_col)
{
    time_step = 0;

    src_row = new_src_position_row + pml_width;
    src_col = new_src_position_col + pml_width;

    // Init field arrays.
    std::fill_n(&dz[0][0], rows * cols, 0.0);
    std::fill_n(&ez[0][0], rows * cols, 0.0);
    std::fill_n(&hx[0][0], rows * cols, 0.0);
    std::fill_n(&hy[0][0], rows * cols, 0.0);
    std::fill_n(&gaz[0][0], rows * cols, 1.0);

    // Init pml arrays.
    std::fill_n(&gi2[0], rows, 1.0);
    std::fill_n(&gi3[0], rows, 1.0);
    std::fill_n(&fi1[0], rows, 0.0);
    std::fill_n(&fi2[0], rows, 1.0);
    std::fill_n(&fi3[0], rows, 1.0);

    std::fill_n(&gj2[0], rows, 1.0);
    std::fill_n(&gj3[0], rows, 1.0);
    std::fill_n(&fj1[0], rows, 0.0);
    std::fill_n(&fj2[0], rows, 1.0);
    std::fill_n(&fj3[0], rows, 1.0);

    // Fill pml arrays.
    for (size_t i = 0; i < pml_width; ++i)
    {
        int xnum = pml_width - i;
        double xd = pml_width;
        double xxn = xnum / xd;
        double xn = 0.33 * std::pow(xxn, 3);

        gi2[i] = 1.0 / (1.0 + xn);
        gi2[rows - 1 - i] = 1.0 / (1.0 + xn);
        gi3[i] = (1.0 - xn) / (1.0 + xn);
        gi3[rows - i - 1] = (1.0 - xn) / (1.0 + xn);

        gj2[i] = 1.0 / (1.0 + xn);
        gj2[rows - 1 - i] = 1.0 / (1.0 + xn);
        gj3[i] = (1.0 - xn) / (1.0 + xn);
        gj3[rows - i - 1] = (1.0 - xn) / (1.0 + xn);

        xxn = (xnum - 0.5) / xd;
        xn = 0.33 * std::pow(xxn, 3);

        fi1[i] = xn;
        fi1[rows - 2 - i] = xn;
        fi2[i] = 1.0 / (1.0 + xn);
        fi2[rows - 2 - i] = 1.0 / (1.0 + xn);
        fi3[i] = (1.0 - xn) / (1.0 - xn);
        fi3[rows - 2 - i] = (1.0 - xn) / (1.0 + xn);

        fj1[i] = xn;
        fj1[rows - 2 - i] = xn;
        fj2[i] = 1.0 / (1.0 + xn);
        fj2[rows - 2 - i] = 1.0 / (1.0 + xn);
        fj3[i] = (1.0 - xn) / (1.0 - xn);
        fj3[rows - 2 - i] = (1.0 - xn) / (1.0 + xn);
    }

    // Fill medium array.
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {

            if (i > pml_width && i < (rows - pml_width) && j > pml_width && j < (cols - pml_width))
            {
                gaz[i][j] = 1.0 / (eps[i - pml_width][j - pml_width] + (sigma[i - pml_width][j - pml_width] * dt) / epsz);
            }
            else
            {
                gaz[i][j] = 1.0 / (epsilon1 + (sigma1 * dt) / epsz);
            }

            // Medium 1.
            // gaz[i][j] = 1.0 / (epsilon1 + (sigma1 * dt) / epsz);

            // Medium 2.
            // if (i>=80 && i<=90 && j>=80 && j<=120) {
            //    gaz[i][j] = 1.0 / (epsilon2 + (sigma2 * dt) / epsz);
            // }

            // if (i>=120 && i<=130 && j>=80 && j<=120) {
            //    gaz[i][j] = 1.0 / (epsilon2 + (sigma2 * dt) / epsz);
            // }

            // // Dielectric border. Medium 0.
            // if (i >= (src_row - 20) && i <= (src_row + 20) && (j == (src_col - 10) || j == (src_col + 10)) || j >= (src_col - 10) && j <= (src_col + 10) && (i == (src_row - 10)))
            // {
            //     gaz[i][j] = 1.0 / (epsilon0 + (sigma0 * dt) / epsz);
            // }
        }
    }
}

// One time layer calculation.
void FdtdPml2D::Calculation()
{

    // Dz field calculation.
    for (int i = 1; i < rows; i++)
    {
        for (int j = 1; j < cols; j++)
        {
            dz[i][j] = gi3[i] * gj3[j] * dz[i][j] + gi2[i] * gj2[j] * 0.5 * (hy[i][j] - hy[i - 1][j] - hx[i][j] + hx[i][j - 1]);
        }
    }

    // Ez field calculation.
    for (int i = 1; i < rows; i++)
    {
        for (int j = 1; j < cols; j++)
        {
            ez[i][j] = gaz[i][j] * dz[i][j];
        }
    }

    // Put Gaussian beam source.
    double source = -2.0 * ((time_step - t0) / tau) * std::exp(-1.0 * std::pow((time_step - t0) / tau, 2));
    ez[src_row][src_col] = source;

    for (int i = 0; i < rows - 1; i++)
    {
        for (int j = 0; j < cols - 1; j++)
        {

            // Hx field calculation.
            hx[i][j] = fj3[j] * hx[i][j] + fj2[j] * 0.5 * (ez[i][j] - ez[i][j + 1]);
            // Hy field calculation.
            hy[i][j] = fi3[i] * hy[i][j] + fi2[i] * 0.5 * (ez[i + 1][j] - ez[i][j]);
        }
    }
}

size_t FdtdPml2D::GetStep()
{
    return 1;
}

size_t FdtdPml2D::GetCurrentTimeStep()
{
    return time_step;
}

FdtdPml2D::FdtdPml2D(std::vector<std::vector<double>> &eps,
                     std::vector<std::vector<double>> &mu,
                     std::vector<std::vector<double>> &sigma,
                     size_t new_src_position_row,
                     size_t new_src_position_col)
{
    SetParams(eps, mu, sigma, new_src_position_row, new_src_position_col);
}

void FdtdPml2D::CalcNextLayer(std::vector<double> &vectX,
                              std::vector<double> &vectY,
                              std::vector<double> &vectEz,
                              std::vector<double> &vectHy,
                              std::vector<double> &vectHx,
                              std::vector<double> &vectEnergy,
                              double &max, double &min)
{

    Calculation();

    size_t step = GetStep();

    // for (int xx = 1; xx < rows - 1; xx += step)
    for (int xx = pml_width; xx < rows - pml_width; xx += step)
    {
        // vectX.push_back(xx);
        // vectY.push_back(xx);
        // for (int yy = 1; yy < cols - 1; yy += step)
        for (int yy = pml_width; yy < cols - pml_width; yy += step)
        {
            // Energy
            // double energy = yy1[xx][yy] * yy1[xx][yy] * Ez1[xx][yy] * Ez1[xx][yy] +
            // Hy1[xx][yy] * Hy1[xx][yy] + Hx1[xx][yy] * Hx1[xx][yy];
            double energy = 1;

            vectX.push_back(xx - pml_width);
            vectY.push_back(yy - pml_width);
            vectEz.push_back(ez[yy][xx]);
            vectHy.push_back(hy[yy][xx]);
            vectHx.push_back(hx[yy][xx]);
            vectEnergy.push_back(energy);

            if (xx > (pml_width + 10) && xx < (rows - pml_width - 10) && yy > (pml_width + 10) && yy < (rows - pml_width - 10))
            {
                if (ez[yy][xx] > max)
                {
                    max = ez[yy][xx];
                }
                if (ez[yy][xx] < min)
                {
                    min = ez[yy][xx];
                }
            }
        }
    }

    time_step++;
}