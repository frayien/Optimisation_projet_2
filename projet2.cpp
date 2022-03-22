#include <iostream>
#include <array>
#include <functional>
#include <optional>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <string>

std::array<std::function<float(float, float)>, 5> P_turb =
{
    [](float x, float y) -> float
    {
        const float p00 =      -233.3;
        const float p10 =       13.44;
        const float p01 =      0.1889;
        const float p20 =     -0.1935;
        const float p11 =    -0.02236;
        const float p02 =    0.005538;
        const float p21 =   0.0004944;
        const float p12 =  -3.527e-05;
        const float p03 =  -1.594e-05;
        return p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p21*x*x*y + p12*x*y*y + p03*y*y*y;
    },
    [](float x, float y) -> float
    {
        const float p00 =      -1.382;
        const float p10 =     0.09969;
        const float p01 =      -1.945;
        const float p20 =   -0.001724;
        const float p11 =     0.09224;
        const float p02 =    0.007721;
        const float p21 =   -0.001096;
        const float p12 =  -6.622e-05;
        const float p03 =  -1.933e-05;
        return p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p21*x*x*y + p12*x*y*y + p03*y*y*y;
    },
    [](float x, float y) -> float
    {
        const float p00 =      -102.4;
        const float p10 =       5.921;
        const float p01 =     -0.5012;
        const float p20 =    -0.08557;
        const float p11 =     0.02467;
        const float p02 =    0.003842;
        const float p21 =  -0.0002079;
        const float p12 =  -2.209e-05;
        const float p03 =  -1.179e-05;
        return p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p21*x*x*y + p12*x*y*y + p03*y*y*y;
    },
    [](float x, float y) -> float
    {
        const float p00 =      -514.8;
        const float p10 =       29.72;
        const float p01 =       2.096;
        const float p20 =     -0.4288;
        const float p11 =     -0.1336;
        const float p02 =    0.005654;
        const float p21 =    0.002048;
        const float p12 =  -5.026e-07;
        const float p03 =  -1.999e-05;
        return p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p21*x*x*y + p12*x*y*y + p03*y*y*y;
    },
    [](float x, float y) -> float
    {
        const float p00 =      -212.1;
        const float p10 =       12.17;
        const float p01 =    0.004397;
        const float p20 =     -0.1746;
        const float p11 =   -0.006808;
        const float p02 =    0.004529;
        const float p21 =   0.0002936;
        const float p12 =  -4.211e-05;
        const float p03 =  -1.176e-05;
        return p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p21*x*x*y + p12*x*y*y + p03*y*y*y;
    },
};

struct params
{
    int Qmax;
    float Hchute;
    std::array<bool, 5> turbines_disponibles;
};

std::ostream & operator<<(std::ostream & out, const params & par)
{
    out << "Qmax : "   << par.Qmax   << std::endl;
    out << "disponibles : [" 
        << par.turbines_disponibles[0] << ", "
        << par.turbines_disponibles[1] << ", "
        << par.turbines_disponibles[2] << ", "
        << par.turbines_disponibles[3] << ", "
        << par.turbines_disponibles[4] << "]"
        << std::endl;
    out << "Hchute : " << par.Hchute << std::endl;
    return out;
}

struct result
{
    std::array<int, 5> Q;
    std::array<float, 5> P;

    int   Qtot() const { return Q[0] + Q[1] + Q[2] + Q[3] + Q[4]; }
    float Ptot() const { return P[0] + P[1] + P[2] + P[3] + P[4]; }
};

std::ostream & operator<<(std::ostream & out, const result & res)
{
    out << "Q : [" 
        << res.Q[0] << ", "
        << res.Q[1] << ", "
        << res.Q[2] << ", "
        << res.Q[3] << ", "
        << res.Q[4] << "] tot : "
        << res.Qtot() << " m^3/s"
        << std::endl;
    out << "P : ["
        << res.P[0] << ", "
        << res.P[1] << ", "
        << res.P[2] << ", "
        << res.P[3] << ", "
        << res.P[4] << "] tot : "
        << res.Ptot() << " MW"
        << std::endl;

    return out;
}

result optimise(params par)
{
    std::unordered_map<int, result> res_old;
    std::unordered_map<int, result> res_new;

    // Pour chaque turbine
    for(std::size_t i = 0; i < 5; ++i)
    {
        // si la turbine est active
        if(!par.turbines_disponibles[i]) continue;
        
        // pour chaque débit Q possible
        for(int Q = 0; Q <= par.Qmax; Q += 5)
        {
            // trouver quel est le débit Qi optimal à faire passer dans la turbine i
            // (le reste du débit, Q - Qi, passera donc dans les turbines d'id < i) 
            result res_max{};
            for(int Qi = 0; Qi <= Q; Qi += 5)
            {
                result res_tmp = res_old[Q - Qi];

                // ajouter la turbine courrante à la solution
                res_tmp.Q[i] = Qi;
                res_tmp.P[i] = P_turb[i](par.Hchute, Qi);

                if(res_tmp.Ptot() > res_max.Ptot())
                {
                    res_max = res_tmp;
                }
            }

            res_new[Q] = res_max;
        }

        res_old = res_new;
    }

    return res_new[par.Qmax];
}

void sanitize(std::string & str)
{
    const std::size_t pos = str.find_first_of(',');
    if(pos != std::string::npos) str.replace(pos, 1, ".");
}

int round_q(float Q)
{
    return static_cast<int>(std::ceil(Q/5.f)) * 5;
}

struct data
{
    int Qmax;
    float Hchute;
    std::array<int, 5> Q;
    std::array<float, 5> P;

    int   Qtot() const { return Q[0] + Q[1] + Q[2] + Q[3] + Q[4]; }
    float Ptot() const { return P[0] + P[1] + P[2] + P[3] + P[4]; }
};

std::ostream & operator<<(std::ostream & out, const data & dat)
{
    out << "Qmax : " << dat.Qmax << std::endl;
    out << "Hchute : " << dat.Hchute << std::endl;
    out << "Q : [" 
        << dat.Q[0] << ", "
        << dat.Q[1] << ", "
        << dat.Q[2] << ", "
        << dat.Q[3] << ", "
        << dat.Q[4] << "] tot : "
        << dat.Qtot() << " m^3/s"
        << std::endl;
    out << "P : ["
        << dat.P[0] << ", "
        << dat.P[1] << ", "
        << dat.P[2] << ", "
        << dat.P[3] << ", "
        << dat.P[4] << "] tot : "
        << dat.Ptot() << " MW"
        << std::endl;

    return out;
}

namespace fs = std::filesystem;

std::vector<data> read_csv(fs::path path, std::size_t to_read)
{
    std::vector<data> lines{};
    lines.reserve(to_read);

    if(!fs::is_regular_file(path)) return lines;

    std::ifstream in(path);

    if(!in.is_open()) return lines;

    const std::uintmax_t FILE_SIZE = fs::file_size(path);

    // begin reading //

    // skip 3 header lines
    in.ignore(FILE_SIZE, '\n');
    in.ignore(FILE_SIZE, '\n');
    in.ignore(FILE_SIZE, '\n');
    
    for(std::size_t i = 0; i<to_read; ++i)
    {
        std::string Eval_str;
        std::string Qtot_str;
        std::string Namont_str;
        std::array<std::string, 5> P_str;
        std::array<std::string, 5> Q_str;

        in.ignore(FILE_SIZE, ';'); // date
        std::getline(in, Eval_str, ';'); // Elav
        std::getline(in, Qtot_str, ';'); // Qtot
        in.ignore(FILE_SIZE, ';'); // Qturb
        in.ignore(FILE_SIZE, ';'); // Qvan
        std::getline(in, Namont_str, ';'); // Namont

        for(std::size_t j = 0; j<5; ++j)
        {
            std::getline(in, Q_str[j], ';'); // Qi
            std::getline(in, P_str[j], ';'); // Pi
        }

        in.ignore(FILE_SIZE, '\n');

        sanitize(Eval_str);
        sanitize(Qtot_str);
        sanitize(Namont_str);

        const int Qmax = round_q(std::stof(Qtot_str));
        const float Hchute = std::stof(Namont_str) - std::stof(Eval_str);

        std::array<int, 5> Q;
        std::array<float, 5> P;

        for(std::size_t j = 0; j<5; ++j)
        {
            sanitize(Q_str[j]);
            sanitize(P_str[j]);
            Q[j] = round_q(std::stof(Q_str[j]));
            P[j] = std::stof(P_str[j]);
        }

        lines.push_back(
        {
            .Qmax = Qmax,
            .Hchute = Hchute,
            .Q = Q,
            .P = P,
        });
    }

    return lines;
}

void write_csv_head(std::ostream & out)
{
    out << "Qtot" << ";" 
        << "Hchute" << ";" 
        << "P1" << ";" << "Q1" << ";" 
        << "P2" << ";" << "Q2" << ";" 
        << "P3" << ";" << "Q3" << ";" 
        << "P4" << ";" << "Q4" << ";" 
        << "P5" << ";" << "Q5" << ";" 
        << ";"
        << "Pr1" << ";" << "Qr1" << ";" 
        << "Pr2" << ";" << "Qr2" << ";" 
        << "Pr3" << ";" << "Qr3" << ";" 
        << "Pr4" << ";" << "Qr4" << ";" 
        << "Pr5" << ";" << "Qr5" << std::endl;
}
void write_csv_line(std::ostream & out, const data & dat, const result & res)
{
    out << dat.Qmax << ";" 
        << dat.Hchute << ";" 
        << res.P[1-1] << ";" << res.Q[1-1] << ";" 
        << res.P[2-1] << ";" << res.Q[2-1] << ";" 
        << res.P[3-1] << ";" << res.Q[3-1] << ";" 
        << res.P[4-1] << ";" << res.Q[4-1] << ";" 
        << res.P[5-1] << ";" << res.Q[5-1] << ";" 
        << ";"
        << dat.P[1-1] << ";" << dat.Q[1-1] << ";" 
        << dat.P[2-1] << ";" << dat.Q[2-1] << ";" 
        << dat.P[3-1] << ";" << dat.Q[3-1] << ";" 
        << dat.P[4-1] << ";" << dat.Q[4-1] << ";" 
        << dat.P[5-1] << ";" << dat.Q[5-1] << std::endl;
}

int main()
{
    fs::path path_in("./DataProjet2022_v1.csv");
    std::vector<data> file_data = read_csv(path_in, 100);

    fs::path path_out("data_out.csv");
    for(int i = 1; fs::exists(path_out); ++i)
    {
        path_out.replace_filename("data_out(" + std::to_string(i) + ").csv");
    }

    std::ofstream out(path_out, std::ios::trunc);

    write_csv_head(out);

    for(const data & line : file_data)
    {
        params par =
        {
            .Qmax = line.Qmax,
            .Hchute = line.Hchute,
            .turbines_disponibles = {true, true, true, true, true}, 
        };

        result res = optimise(par);

        write_csv_line(out, line, res);
    }

    out.close();

    return 0;
}