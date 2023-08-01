#include <bits/stdc++.h>
#include <fstream>
// #include "gnuplot-iostream.h"
// #include "gnuplot.h"

#include <chrono>
using namespace std::chrono;

using namespace std;
using ll = long long;
#define lb lower_bound
#define ub upper_bound
#define all(v) (v).begin(), (v).end()
#define all1(v) (v).begin() + 1, (v).end()
#define allr(v) (v).rbegin(), (v).rend()
#define allr1(v) (v).rbegin() + 1, (v).rend()
#define sort0(v) sort(all(v))
#define fo(i, a, b) for (i = a; i <= b; i++)
#define fi(i, a, b) for (i = a; i >= b; i--)

// freopen("input.txt", "r", stdin);
// freopen("output.txt", "w", stdout);
#define rep(i, begin, end) for (__typeof(end) i = (begin) - ((begin) > (end)); i != (end) - ((begin) > (end)); i += 1 - 2 * ((begin) > (end)))
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef pair<ll, ll> pll;

#define sz(x) (ll) x.size()
#define pb push_back
#define ppb pop_back
#define mkp make_pair
#define inf 1000000000000000005

const ll mod = 1e9 + 7;
const ll mod1 = 998244353;
const ll N = 1e5; // limit for array size

class Plate_Fin_variables
{
public:
    double V, A, Tf, Ab, Tin, Af, Ap, As, Auf, Atot, P, L1, W1, B, H, t, L, W, Tb, Tinf, N, S, qb, Tmax;
    Plate_Fin_variables()
    {
        Tinf = 298;
        Tb = 373;
        Tin = 350;
        Tmax = 400;
        V = 1;
        H = 900 * 0.001;
        L = 200 * 0.001;
        L1 = 1;
        W1 = 1;
        W = 47 * 0.001;
        t = 10 * 0.001;
        N = 2;
        S = 4 * 0.001;
        qb = 20000;

        solve();
    }

    void change_fin_dim(double L, double H, double N, double t)
    {
        this->L = L;
        this->H = H;
        // this->W = W;
        this->t = t;
        this->N = N;
        // this->S = S;
        solve();
    }
    void solve()
    {
        P = 2 * (t + L);
        A = W1 * L;
        Af = L * t;
        As = P * H;
        Ap = t * H;
        W1 = V / (L * H);
        S = (W1 - (N * t)) / (N - 1);
        Auf = L * W1 - (N * Af);
        Atot = Auf + As;
        Tf = (Tinf + Tb) / 2;
    }
};
class Pin_Fin_variables
{
public:
    double A, Tf, Ab, Af, Ap, As, Auf, Atot, P, L1, W1, B, H, D, L, W, Tb, Tinf, N, S, qb;
    Pin_Fin_variables()
    {
        Tinf = 293;
        Tb = 373;
        H = 100 * 0.001;
        L = 100 * 0.001;
        L1 = 1;
        W1 = 1;
        W = 100 * 0.001;
        D = 10 * 0.001;
        N = 625;
        S = 4;
        qb = 10;

        solve();
    }

    void change_fin_dim(double H, double L, double D, double W, double N, double S)
    {
        this->L = L;
        this->H = H;
        this->D = D;
        this->W = W;
        this->N = N;
        this->S = S;
        solve();
    }
    void solve()
    {
        P = M_PI * D;
        A = W1 * L1;
        Af = M_PI * D * D / 4;
        As = P * H;
        // Ap = t * H;
        Auf = L * W - (N * Af);
        Atot = Auf + As;
        Tf = (Tinf + Tb) / 2;
    }
};

class Flow_Parameters
{
public:
    double k, ka, NuL, Nux, Rex, Pe, Pr, q, RaL, ReL, mu, nu, h, qs, uinf, rho, hx, havg;
    Flow_Parameters()
    {

        qs = 0;
        k = 247;
        h = 80;
        uinf = 1.6;
        hx = 80;
        solve();
    }

    void solve()
    {
        Pr = 0.7202;
        rho = 1.06;
        nu = 1.896 * 1e-5;
        mu = 2.008 * 1e-5;
        ka = 0.02808;

        // Nux = 0.453
    }
    void solve_plateh(double x)
    {
        Rex = uinf * x / nu;
        if (Rex < (5 * 1e5))
            Nux = 0.332 * pow(Rex, (double)1 / 2) * pow(Pr, (double)1 / 3);
        else
            Nux = 0.0296 * pow(Rex, (double)0.8) * pow(Pr, (double)1 / 3);

        hx = Nux * ka / x;
    }
    void solve_pinh(double D, double x)
    {
        Rex = uinf * D / nu;
        // Nux = 0.023 * pow(Rex, (double)0.8) * pow(Pr, (double)1 / 3);
        Nux = 2 + (0.6 * pow(Rex, 0.5) * pow(Pr, (double)1 / 3) * pow((D / x), 0.5)) / pow((1 + 0.4 * (pow(Rex, 0.5) * pow(Pr, (double)1 / 3) * pow((D / x), 0.5))), 0.36);
        // Nux = 1.62 * pow(Rex, 0.54) * pow((y / D), -1.53) * pow((x / D), -0.48) * pow((H / D), 0.78);
        hx = Nux * ka / D;
    }
};
vector<double> TDMA1D(ll n, vector<vector<double>> &a)
{
    ll i, j, m;
    m = 4;
    vector<double> phi(n + 2, 0);

    vector<double> P(n + 1, 0), Q(n + 1, 0);
    fo(i, 1, n)
    {
        P[i] = a[i][2] / (a[i][1] - a[i][3] * P[i - 1]);
        Q[i] = (a[i][3] * Q[i - 1] + a[i][4]) / (a[i][1] - a[i][3] * P[i - 1]);
    }
    fi(i, n, 1)
    {
        phi[i] = P[i] * phi[i + 1] + Q[i];
    }
    return phi;
}

vector<vector<double>> TDMA2D(ll n, ll m, vector<vector<vector<double>>> &a, ll sweep)
{
    ll i, j, k;
    double err = 1e9;
    ll sw = sweep;
    if (sweep == 3 || sweep == 5)
    {
        sweep--;
    }
    if (sweep > 3)
    {
        // sweep--;
        swap(n, m);
    }
    ll iter = 0;
    vector<vector<double>> phi(n + 2, vector<double>(m + 2, 0));
    while (err >= 1e-2 && iter <= 2000)
    {
        err = 0.0;
        fo(i, 1, m)
        {
            ll ind = (sw == 2 || sw == 4) ? i : m - i + 1;
            vector<vector<double>> tda(n + 1, vector<double>(5, 0));

            for (j = 1; j <= n; j++)
            {
                double val = 0;
                ll f1 = j, f2 = ind;
                if (sw > 3)
                {
                    f1 = ind;
                    f2 = j;
                }
                for (k = 1; k < 6; k++)
                {
                    if (k == sweep)
                    {
                        val += a[f1][f2][k] * ((sw <= 3) ? phi[f1][f2 + 1] : phi[f1 + 1][f2]);
                        continue;
                    }
                    else if (k == sweep + 1)
                    {
                        val += a[f1][f2][k] * ((sw <= 3) ? phi[f1][f2 - 1] : phi[f1 - 1][f2]);
                        continue;
                    }
                    tda[j][max(1LL, k - ((sw < 4) ? 2 : 0))] = a[f1][f2][k];
                }
                tda[j][4] = val + a[f1][f2][6];
            }

            vector<double> ph1 = TDMA1D(n, tda);

            fo(j, 1, n)
            {
                ll f1 = j, f2 = ind;
                if (sw > 3)
                {
                    f1 = ind;
                    f2 = j;
                }
                err += (phi[f1][f2] - ph1[j]) * (phi[f1][f2] - ph1[j]);
                phi[f1][f2] = ph1[j];
            }
        }
        err = sqrt(err);

        iter++;
    }
    if (iter > 1e5)
    {
        cout << "Diverges!!!\n";
    }
    // for (i = 1; i <= m; i++)
    // {
    //     for (j = 1; j <= n; j++)
    //     {
    //         cout << setprecision(4) << fixed << phi[j][i] << " ";
    //     }
    //     cout << "\n";
    // }
    return phi;
}

void fin_1D_plate_sample()
{
    ofstream PCT("plate_central_temp.txt");
    ofstream PH("plate_heat.txt");
    ofstream PR("plate_Resistance.txt");
    ofstream NOF("plate_nof.txt");
    // ofstream BT("plate_base_temp.txt");
    // Gnuplot gp;

    ll n, i, m, j, cl;

    // Insulated B.C
    Plate_Fin_variables f1;
    Flow_Parameters c1;
    cin >> n >> m;
    ll change_H, change_N, change_t;
    double del_H = 50 * 0.001, del_N = 18, O_N = f1.N, O_H = f1.H, H_1 = f1.H, O_t = f1.t, del_t = 1 * 0.001;
    fo(cl, 1, 50)
    {
        f1.change_fin_dim(f1.L + (cl > 0 ? (100 * 0.001) : 0), O_H, O_N, O_t);
        PH << f1.L << "\n[";
        PR << f1.L << "\n[";
        // BT << f1.L << "\n";
        NOF << f1.L << "\n";

        fo(change_H, 0, 1)
        {
            H_1 = f1.H + (change_H > 0 ? 1 : 0) * del_H;
            f1.change_fin_dim(f1.L, H_1, O_N, O_t);
            PH << "[";
            PR << "[";

            fo(change_N, 1, 1)
            {

                double N_1 = f1.N + (change_N > 0 ? 1 : 0) * del_N;
                f1.change_fin_dim(f1.L, f1.H, N_1, O_t);
                fo(change_t, 0, 1)
                {
                    double t_1 = f1.t + (change_t > 0 ? 1 : 0) * del_t;
                    f1.change_fin_dim(f1.L, f1.H, f1.N, t_1);
                    vector<vector<vector<double>>> v(n + 1, vector<vector<double>>(m + 1, vector<double>(7)));
                    double del_x = f1.L / n, del_y = f1.H / m;

                    fo(i, 1, n)
                    {

                        double x_conv = del_x * (i - 0.5);
                        c1.solve_plateh(x_conv);

                        // cout << c1.hx << "\n";
                        fo(j, 1, m)
                        {
                            double ap = 0, aw = del_y * c1.k * f1.t / del_x, ae = del_y * c1.k * f1.t / del_x, an = del_x * c1.k * f1.t / del_y, as = del_x * c1.k * f1.t / del_y, ab = 0, b = 0;

                            double P = 2 * del_y, Ueq = c1.hx * (2 * c1.k / del_x) / (c1.hx + (2 * c1.k / del_x)), msq = c1.hx * 2 * (del_x + del_y) / (del_x * del_y);
                            if (i == 1)
                            {
                                ab += 2 * as;
                                as = 0;
                                b += ab * f1.Tb;
                                // b += f1.qb * del_x * f1.t;
                            }
                            else if (i == n)
                            {
                                an = 0;
                            }

                            if (j == 1)
                            {
                                ab += Ueq * f1.t * del_y;
                                aw = 0;
                                b += Ueq * f1.Tinf * f1.t * del_y;
                            }
                            else if (j == m)
                            {
                                ab += Ueq * f1.t * del_y;
                                ae = 0;
                                b += Ueq * f1.Tinf * f1.t * del_y;
                            }

                            ap = aw + ae + an + as + ab + (msq * del_x * del_y * f1.t);
                            b += msq * del_x * del_y * f1.t * f1.Tinf;
                            v[i][j][1] = ap, v[i][j][2] = ae, v[i][j][3] = aw, v[i][j][4] = an, v[i][j][5] = as, v[i][j][6] = b;
                        }
                    }
                    vector<vector<double>> phi = TDMA2D(n, m, v, 4);

                    // Temperature and heat flux at different points

                    double Qf = 0, Qmf = 0, Qt = 0, Quf = 0, havg = 0;
                    for (i = 1; i <= n; i++)
                    {
                        double x_conv = del_x * (i - 0.5);
                        c1.solve_plateh(x_conv);
                        havg += c1.hx;
                        double Tinnn = 0;
                        PCT << "[";
                        for (j = 1; j <= m; j++)
                        {
                            double y = i * del_x, x = j * del_y;
                            double P = 2 * (del_x);
                            if (j == 1 || j == m)
                                P += f1.t;
                            PCT << setprecision(4) << fixed << phi[i][j] << ",";

                            Qf += c1.hx * P * del_y * (phi[i][j] - f1.Tinf);
                        }
                        PCT << "],\n";
                    }
                    PCT << "\n";
                    // cout << Qf << "\n";
                    for (i = 1; i <= n; i++)
                    {
                        double x_conv = del_x * (i - 0.5);
                        c1.solve_plateh(x_conv);
                        Qmf += c1.hx * f1.P * del_y * (f1.Tb - f1.Tinf);
                    }
                    havg /= n;
                    // cout << havg << "\n";
                    c1.havg = havg;
                    Quf += havg * f1.Auf * (f1.Tb - f1.Tinf);

                    double nf = Qf / Qmf; // c1.hx * f1.Auf
                    Qt = (nf * f1.N * Qf) + (Quf);
                    double Rth = (f1.Tb - f1.Tinf) / Qt;
                    // cout << Qf << " " << nf << " " << Qt << " " << Rth << "\n";
                    if (f1.S > 0)
                    {
                        PH << setprecision(4) << fixed << Qt << ",";
                        // BT << setprecision(4) << fixed << f1.Tb << ",";
                        PR << setprecision(10) << fixed << Rth << ",";
                    }
                    else
                    {
                        PH << setprecision(4) << fixed << f1.S << " -1"
                           << "\n";
                        // BT << setprecision(4) << fixed << f1.Tb << ",";
                        PR << setprecision(4) << fixed << f1.S << " -1"
                           << "\n";
                        NOF << setprecision(4) << fixed << f1.S << " -1"
                            << "\n";
                    }
                }
                // PH << ",";
                // // BT << "\n";
                // PR << ",";
                NOF << setprecision(4) << fixed << f1.N << ",";
            }
            PH << "]\n";
            // BT << "\n";
            PR << "]\n";
            NOF << "\n";
        }
        PH << "]\n";
        // BT << "\n";
        PR << "]\n";
        NOF << "\n";
    }

    PCT.close();
    // BT.close();
    PH.close();
    PR.close();
    NOF.close();
}

// void test_plate()
// {
//     ll i;
//     double m = sqrt((c1.havg * f1.P) / (c1.k * f1.Af));
//     double del_y = f1.H / 50, y = 0;
//     fo(i, 1, 50)
//     {
//         y += del_y;
//         double T = f1.Tinf + (f1.Tb - f1.Tinf) * (cosh(m * (f1.H - y)) / cosh(m * (f1.H)));
//         cout << T << ", ";
//     }
//     cout << "\n\n";
//     double Qf = m * (c1.k * f1.Af) * (f1.Tb - f1.Tinf) * tanh(m * f1.H), eta = Qf / (c1.havg * f1.As * (f1.Tb - f1.Tinf)), Quf = c1.havg * (f1.Auf * (f1.Tb - f1.Tinf));
//     double Qt = (Qf * eta * f1.N + Quf), Rth = (f1.Tb - f1.Tinf) / (Qt);
//     cout << Qf << " " << eta << " " << Qt << " " << Rth << "\n";
// }
void fin_1D_pin_sample()
{
    ofstream MyFile("pin_fin.txt");

    ll n, i, m, j, row;
    Pin_Fin_variables f1;
    Flow_Parameters c1;
    // Insulated B.C
    cin >> n >> m;
    ll rows = sqrt(f1.N);
    vector<vector<vector<double>>> v(n + 1, vector<vector<double>>(m + 1, vector<double>(7)));
    vector<vector<double>> temp, phi(rows + 1, vector<double>(n + 1));
    double del_z = f1.H / n, del_r = f1.D / (2 * m), del_dis = f1.L / (rows + 1);
    fo(row, 1, rows)
    {

        double x_conv = del_dis * (row);
        c1.solve_pinh(f1.D, x_conv);
        cout << c1.hx << "\n";

        fo(i, 1, n)
        {
            fo(j, 1, m)
            {
                double rn = del_r * (i - 0.5), rs = del_r * (i - 0.5), rw = del_r * (i - 1), re = del_r * (i);
                double ap = 0, aw = rw * c1.k * del_z / (del_r), ae = re * c1.k * del_z / (del_r), an = rn * del_r * c1.k / del_z, as = rs * del_r * c1.k / del_z, ab = 0, b = 0;

                double Ueq = c1.hx * (2 * c1.k / del_r) / (c1.hx + (2 * c1.k / del_r));
                if (i == 1)
                {

                    ab += 2 * as;
                    as = 0;
                    b += ab * f1.Tb;
                    // b += f1.qb * del_x * f1.t;
                }
                else if (i == n)
                {
                    an = 0;
                }

                if (j == 1)
                {
                    aw = 0;
                }
                else if (j == m)
                {
                    re = f1.D / 2;
                    ab += Ueq * re * del_z;
                    ae = 0;
                    b += Ueq * f1.Tinf * re * del_z;
                }

                ap = aw + ae + an + as + ab;
                v[i][j][1] = ap, v[i][j][2] = ae, v[i][j][3] = aw, v[i][j][4] = an, v[i][j][5] = as, v[i][j][6] = b;
            }
        }
        temp = TDMA2D(n, m, v, 4);
        fo(i, 1, n)
        {
            phi[row][i] = temp[i][m];
        }
    }

    // Temperature and heat flux at different points
    double Qf = 0, Qmf = 0, Qt = 0, Quf = 0, havg = 0;
    for (i = 1; i <= rows; i++)
    {

        double x_conv = del_dis * (i);
        c1.solve_pinh(f1.D, x_conv);
        Qmf += c1.hx * f1.As * (f1.Tb - f1.Tinf);

        havg += c1.hx;
        for (j = 1; j <= n; j++)
        {
            double P = f1.P;

            MyFile << setprecision(4) << fixed << phi[i][j] << " " << c1.hx * P * del_z * (phi[i][j] - f1.Tinf) << "  ";
            Qf += c1.hx * P * del_z * (phi[i][j] - f1.Tinf);
        }
        MyFile << "\n";
    }
    havg /= rows;
    Quf = havg * f1.Auf * (f1.Tb - f1.Tinf);

    double nf = Qf / Qmf; // c1.hx * f1.Auf
    Qt = (nf * rows * Qf) + (Quf);
    // cout << Qmf << " " << Quf << "  ";
    cout << Qf << " " << nf << " " << Qt << "\n";
}
int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    auto start = high_resolution_clock::now();

    fin_1D_plate_sample();
    // test_plate();
    // fin_1D_pin_sample();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "\n\n"
         << duration.count() * 1e-6 << " seconds" << endl;

    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)
