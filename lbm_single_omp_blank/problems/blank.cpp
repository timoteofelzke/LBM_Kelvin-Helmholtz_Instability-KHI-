/**
 *  \file default.cpp
 *
 *  Description: A empty problem (default) to be used as template to implement new problems
 *  into the code
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "../lbmbase/geometry.h"
#include "../parameters.h"
#include "../lbmbase/boundary/neumann.h"
#include "../lbmbase/lattice/lattice.h"

#include "blank.h"

using namespace std;

namespace blankProblem
{
/*
stProblem problemParameters( string filename )
{
    stProblem prm;
    using boost::property_tree::ptree;
    ptree pt;
    read_ini(filename.c_str(),pt);
    return prm;
}
*/

/*stProblem problemParameters(string filename) {

    stProblem prm;
    try {
        using boost::property_tree::ptree;
        ptree pt;
        
        // Verificação básica do arquivo
        ifstream file(filename);
        if (!file) {
            throw runtime_error("Arquivo de configuração não encontrado: " + filename);
        }
        file.close();

        // Leitura do arquivo INI
        read_ini(filename, pt);

        // Atribuição segura de parâmetros com valores padrão
        try {
            prm.tau = pt.get<double>("BGK.tau");  // Valor prioritário do BGK
        } catch (...) {
            prm.tau = pt.get<double>("General.tau", 0.6);  // Fallback para General.tau
        }

        // Parâmetros essenciais com fallback
        prm.initial_density = pt.get<double>("General.initial_density", 1.0);
        prm.steps = pt.get<int>("General.number_of_steps", 1500);
        
        // Geometria (necessária para problemInitialize)
        prm.size_x = pt.get<int>("Geometry.size_x", 256);
        prm.size_y = pt.get<int>("Geometry.size_y", 256);

    } catch (const boost::property_tree::ptree_bad_data& e) {
        cerr << "ERRO: Valor inválido no arquivo " << filename << endl;
        cerr << "Detalhe: " << e.what() << endl;
        cerr << "Verifique:" << endl;
        cerr << "1. Uso de '.' como separador decimal" << endl;
        cerr << "2. Valores numéricos não vazios" << endl;
        exit(EXIT_FAILURE);
        
    } catch (const exception& e) {
        cerr << "Erro ao processar " << filename << ": " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    
    return prm;
}*/
stProblem problemParameters(string filename) {
    stProblem prm;
    try {
        using boost::property_tree::ptree;
        ptree pt;

        // 1. Verificação rigorosa do arquivo
        ifstream file(filename);
        if (!file.is_open()) {
            throw runtime_error("Arquivo não encontrado: " + filename);
        }

        // 2. Leitura segura do arquivo
        try {
            read_ini(filename, pt);
        } catch (const boost::property_tree::ini_parser_error& e) {
            throw runtime_error("Erro de sintaxe no arquivo: " + string(e.what()));
        }

        // 3. Função auxiliar para conversão segura
        auto safe_get = [&](const string& path, auto default_value) {
            try {
                return pt.get<decltype(default_value)>(path);
            } catch (...) {
                cerr << "Aviso: Usando valor padrão para " << path << endl;
                return default_value;
            }
        };

        // 4. Carregamento dos parâmetros com verificações
        prm.tau = safe_get("BGK.tau", safe_get("General.tau", 0.6));
        prm.initial_density = safe_get("General.initial_density", 1.0);
        prm.steps = safe_get("General.number_of_steps", 1500);
        prm.size_x = safe_get("Geometry.size_x", 256);
        prm.size_y = safe_get("Geometry.size_y", 256);

        // 5. Verificação de valores críticos
        if (prm.tau <= 0.5 || prm.tau >= 2.0) {
            throw runtime_error("Valor de tau inválido. Deve estar entre 0.5 e 2.0");
        }

    } catch (const exception& e) {
        cerr << "ERRO FATAL: " << e.what() << endl;
        cerr << "Verifique o arquivo " << filename << endl;
        cerr << "Certifique-se que todos os valores numéricos:" << endl;
        cerr << "1. Usam ponto (.) como separador decimal" << endl;
        cerr << "2. Não contêm caracteres inválidos" << endl;
        cerr << "3. Estão nas seções corretas" << endl;
        exit(EXIT_FAILURE);
    }

    return prm;
}

/*void problemInitialize(double* iniN, stCollision& col, stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force)
{
    double U0 = 0.0625;                // Velocidade base
    double H = 1.0;                // Controle da espessura da camada (pode ser ajustado) 10.0
    double yc = geo.ny / 2.0;         // Centro da camada em y
    int n_modes = 4; // Modo de perturbação (ajustar conforme necessário)
    double k = 2.0 * M_PI * n_modes / static_cast<double>(geo.nx); //2.0 * M_PI / static_cast<double>(geo.nx); // Número de onda para perturbação k = 2.0 * M_PI / (geo.nx / 2.0);  
    double epsilon = 0.05;          // Amplitude da perturbação 0.01
    double sigma = H/3.0;               // Largura da perturbação
    double rho = 1.0;               // Densidade
    for (int y = 0; y < geo.ny; ++y) {
        for (int x = 0; x < geo.nx; ++x) {
            int idx = (y * geo.nx + x) * NUM_OF_VEL;
            
            // Perfil de cisalhamento suave: U(y)=U0 tanh((y-yc)/H)
            double u = U0 * tanh((y - yc) / H);
            
            // Perturbação vertical: v'(x,y)=ε U0 sin(k x) exp(-((y-yc)^2)/σ^2)
            double v_perturb = epsilon * U0 * sin(k * x) * exp(-((y - yc) * (y - yc)) / (sigma * sigma));
            
            // Considerando que o fluxo primário está na direção horizontal (u)
            // E que a perturbação é aplicada em v
            equilibriumDistribution(u, v_perturb, 0.0, rho, iniN + idx);
        }
    }
}
    */
   void problemInitialize(double* iniN, stCollision& col, stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force)
   {
       double U0 = 0.0625;               // Igual ao Python
       double lbd = 80.0;                // Novo: Controle de espessura (λ do Python)
       double eps = 0.05;                // Amplitude igual ao Python
       double Re = 6000.0;               // Reynolds igual ao Python
       double rho = 1.0;                 
       int L = geo.nx;                   // L = 96 como no Python
   
       // ========== Cálculo de τ idêntico ao Python ==========
       double nu = (U0 * L) / Re;        // ν = 0.0625*96/6000 = 0.001
       col.tau = 3.0 * nu + 0.5;         // τ = 3*0.001 + 0.5 = 0.503 (Python: tau=0.503)
       col.kinematicViscosity = nu;      // Garante compatibilidade
   
       // ========== Perfil de Velocidade com Duas Camadas ==========
       double y_quarter = geo.ny / 4.0;  // Y = NY/4 (primeira camada)
       double y_3quarter = 3.0 * geo.ny / 4.0; // Y = 3NY/4 (segunda camada)
   
       for (int y = 0; y < geo.ny; ++y) {
           for (int x = 0; x < geo.nx; ++x) {
               int idx = (y * geo.nx + x) * NUM_OF_VEL;
               
               // 1. Perfil de velocidade horizontal (igual ao Python)
               double u;
               if (y <= geo.ny/2) {
                   u = U0 * tanh(lbd * (y/(double)geo.ny - 0.25));  // Camada inferior
               } else {
                   u = U0 * tanh(lbd * (0.75 - y/(double)geo.ny));  // Camada superior
               }
               
               // 2. Perturbação vertical idêntica ao Python
               double phase = (x/(double)geo.nx) + 0.25;            // X/NX + 1/4
               double v_perturb = eps * U0 * sin(2.0 * M_PI * phase); 
   
               // 3. Inicialização das distribuições
               equilibriumDistribution(u, v_perturb, 0.0, rho, iniN + idx);
           }
       }
   }

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
    // Finalização do problema (pode ser usado para cálculos finais ou logs)
}

/*
void problemBoundary(double* iniN, stCollision& col, stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step)
{
    // Direção x: Condições periódicas (mantém-se)
    for (int y = 0; y < geo.ny; ++y) {
        int idxL = (y * geo.nx) * NUM_OF_VEL;  // Primeira célula da linha
        int idxR = (y * geo.nx + (geo.nx - 1)) * NUM_OF_VEL;  // Última célula da linha
    
        // Copia da borda direita para a esquerda
        memcpy(iniN + idxL, iniN + idxR, NUM_OF_VEL * sizeof(double));
    
        // Copia da borda esquerda para a direita
        memcpy(iniN + idxR, iniN + idxL, NUM_OF_VEL * sizeof(double));
    }
    
    // Direção y: Aplicando condições free-slip
    // Free-slip em y=0
    for (int x = 0; x < geo.nx; ++x) {
        int idx_top = (x) * NUM_OF_VEL;  // Posição y = 0
        int idx_next = (geo.nx + x) * NUM_OF_VEL; // Posição y = 1 (primeira célula interna)

        // Copia as velocidades tangenciais, mantém v = 0
        memcpy(iniN + idx_top, iniN + idx_next, NUM_OF_VEL * sizeof(double));
        *(iniN + idx_top + 1) = 0.0;  // Garante que vy = 0
    }

    // Free-slip em y = ny-1
    for (int x = 0; x < geo.nx; ++x) {
        int idx_bot = ((geo.ny - 1) * geo.nx + x) * NUM_OF_VEL;  // Posição y = ny-1
        int idx_prev = ((geo.ny - 2) * geo.nx + x) * NUM_OF_VEL;   // Posição y = ny-2

        memcpy(iniN + idx_bot, iniN + idx_prev, NUM_OF_VEL * sizeof(double));
        *(iniN + idx_bot + 1) = 0.0;  // Garante que vy = 0
    }
    
}
*/
/*
void problemBoundary(double* iniN, stCollision& col, stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step)
{
    // ========== Condições Periódicas em X ==========
    for (int y = 0; y < geo.ny; ++y) {
        int idxL = (y * geo.nx) * NUM_OF_VEL;          // Primeira coluna (x=0)
        int idxR = (y * geo.nx + (geo.nx - 1)) * NUM_OF_VEL; // Última coluna (x=nx-1)
        
        // Troca bidirecional para periodicidade perfeita
        double temp[NUM_OF_VEL];
        memcpy(temp, iniN + idxL, NUM_OF_VEL * sizeof(double));
        memcpy(iniN + idxL, iniN + idxR, NUM_OF_VEL * sizeof(double));
        memcpy(iniN + idxR, temp, NUM_OF_VEL * sizeof(double));
    }

    // ========== Condições Free-Slip em Y ==========
    for (int x = 0; x < geo.nx; ++x) {
        // Borda inferior (y = 0)
        int idx_bottom = (0 * geo.nx + x) * NUM_OF_VEL;
        int idx_above = (1 * geo.nx + x) * NUM_OF_VEL;
        
        // Espelhamento das distribuições (D2Q9)
        iniN[idx_bottom + 2] = iniN[idx_above + 4];  // f2 = f4 (direção vertical)
        iniN[idx_bottom + 5] = iniN[idx_above + 7];  // f5 = f7 (diagonal inferior-direita)
        iniN[idx_bottom + 6] = iniN[idx_above + 8];  // f6 = f8 (diagonal inferior-esquerda)
        iniN[idx_bottom + 1] = 0;                  // vy = 0

        // Borda superior (y = ny-1)
        int idx_top = ((geo.ny-1) * geo.nx + x) * NUM_OF_VEL;
        int idx_below = ((geo.ny-2) * geo.nx + x) * NUM_OF_VEL;
        
        iniN[idx_top + 4] = iniN[idx_below + 2];     // f4 = f2
        iniN[idx_top + 7] = iniN[idx_below + 5];     // f7 = f5
        iniN[idx_top + 8] = iniN[idx_below + 6];     // f8 = f6
        iniN[idx_top + 1] = 0;                     // vy = 0
    }
}
*/

void problemBoundary(double* iniN, stCollision& col, stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step)
{
    // ========== Condições Periódicas em X ==========
    for (int y = 0; y < geo.ny; ++y) {
        int idxL = (y * geo.nx) * NUM_OF_VEL;          // Primeira coluna (x=0)
        int idxR = (y * geo.nx + (geo.nx - 1)) * NUM_OF_VEL; // Última coluna (x=nx-1)
        
        // Troca bidirecional para periodicidade perfeita
        double temp[NUM_OF_VEL];
        memcpy(temp, iniN + idxL, NUM_OF_VEL * sizeof(double));
        memcpy(iniN + idxL, iniN + idxR, NUM_OF_VEL * sizeof(double));
        memcpy(iniN + idxR, temp, NUM_OF_VEL * sizeof(double));
    }

    // ========== Condições de Paredes Livres (Open Boundaries) em Y ==========
    for (int x = 0; x < geo.nx; ++x) {
        // Borda inferior (y = 0)
        int idx_bottom = (0 * geo.nx + x) * NUM_OF_VEL;
        int idx_above = (1 * geo.nx + x) * NUM_OF_VEL;
        
        // Extrapolação linear para manter dv/dy = 0
        for (int i = 0; i < NUM_OF_VEL; ++i) {
            iniN[idx_bottom + i] = iniN[idx_above + i];
        }

        // Borda superior (y = ny-1)
        int idx_top = ((geo.ny-1) * geo.nx + x) * NUM_OF_VEL;
        int idx_below = ((geo.ny-2) * geo.nx + x) * NUM_OF_VEL;
        
        for (int i = 0; i < NUM_OF_VEL; ++i) {
            iniN[idx_top + i] = iniN[idx_below + i];
        }
    }
}
void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
    // Atualizações pós-colisão, se necessário
}

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
    return true;
}


}

