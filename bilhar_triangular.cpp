#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// Tolerância numérica (10^-6) usada para evitar erros de arredondamento 
// de ponto flutuante ao verificar colisões e interseções.
const double WIR = 1e-6;

// Função para verificar se o ponto está dentro dos limites do segmento
bool ponto_no_segmento(double x, double y, double x1, double y1, double x2, double y2) {
    // Verifica se é colinear com tolerância numérica
    if (fabs((y2 - y1) * (x - x1) - (y - y1) * (x2 - x1)) > WIR) {
        return false;
    }
    return (x >= fmin(x1, x2) - WIR && x <= fmax(x1, x2) + WIR &&
            y >= fmin(y1, y2) - WIR && y <= fmax(y1, y2) + WIR);
}

// ========================================================================
// Calcula a reflexão especular da partícula usando a fórmula vetorial:
// v_final = v_inicial - 2 * (v_inicial · normal) * normal
// Mantém a magnitude da velocidade constante (choque perfeitamente elástico).
// ========================================================================
void calcular_reflexao(double *vx, double *vy, double a_superficie, double b_superficie) {
    // Vetor normal à superfície
    double nx = -a_superficie;
    double ny = 1.0;
    
    // Normalizar o vetor normal
    double norm = sqrt(nx*nx + ny*ny);
    nx /= norm;
    ny /= norm;
    
    // Produto escalar entre velocidade e normal
    double dot = (*vx) * nx + (*vy) * ny;
    
    // Reflexão: v' = v - 2*(v·n)*n
    *vx = *vx - 2 * dot * nx;
    *vy = *vy - 2 * dot * ny;
    
    // Normalizar para manter velocidade constante
    double v_norm = sqrt((*vx)*(*vx) + (*vy)*(*vy));
    *vx = (*vx) / v_norm;
    *vy = (*vy) / v_norm;
}

// Evita falsas colisões: garante que a nova interseção detectada não seja 
// exatamente o mesmo ponto de onde a partícula acabou de refletir.
bool intersecao_valida(double x_intersec, double y_intersec, double x_inicial, double y_inicial) {
    return fabs(x_intersec - x_inicial) > WIR || fabs(y_intersec - y_inicial) > WIR;
}

int main() {
    // ========================================================================
    // [1] PARÂMETRO DO BILHAR (N)
    // Modifique o valor de N abaixo para alterar os lados do triângulo.
    // Lados: AC = N, BC = N + 1, AB = N + 2
    // ========================================================================
    int N = 7;
    double lado_AB = N + 2;  // Lado maior
    double lado_AC = N;      // Lado menor
    double lado_BC = N + 1;  // Lado intermediário
    
    printf("Triângulo com lados: %.0f, %.0f, %.0f\n", lado_AC, lado_BC, lado_AB);
    
    // Coordenadas dos vértices
    double xa, ya, xb, yb, xc, yc;
    xa = 0; ya = 0;
    xb = lado_AB; yb = 0;
    
    // Cálculo das coordenadas do vértice C usando a Lei dos Cossenos.
    // Isso garante que o triângulo feche perfeitamente com os lados informados.
    double cos_A = (lado_AB*lado_AB + lado_AC*lado_AC - lado_BC*lado_BC) / (2 * lado_AB * lado_AC);
    double sin_A = sqrt(1 - cos_A*cos_A);
    xc = lado_AC * cos_A;
    yc = lado_AC * sin_A;
    
    printf("Vértices:\n");
    printf("A: (0, 0)\n");
    printf("B: (%.6f, 0)\n", lado_AB);
    printf("C: (%.6f, %.6f)\n", xc, yc);
    
    // ========================================================================
    // [2] CONDIÇÕES INICIAIS (POSIÇÃO E VELOCIDADE INICIAL)
    // Modifique as variáveis abaixo para alterar o ponto de partida e a 
    // direção (velocidade inicial) do lançamento da partícula.
    // ========================================================================
    double ponto_inicial_x = lado_AB / 2.0; // Posição X inicial (ex: meio da base)
    double ponto_inicial_y = 0.0;           // Posição Y inicial (na borda AB)
    double angulo_inicial = 89.9;           // Ângulo de lançamento em graus
         
    // A velocidade inicial (vx, vy) é calculada a partir do ângulo inicial
    double angulo_rad = angulo_inicial * M_PI / 180.0;
    double vx = cos(angulo_rad); // Componente X da velocidade inicial
    double vy = sin(angulo_rad); // Componente Y da velocidade inicial
    
    double x = ponto_inicial_x;
    double y = ponto_inicial_y;
    double a = vy/vx;
    double b = y - a*x;
    
    // Coeficientes de cada lado
    double a_AB = 0;
    double b_AB = 0;
    double a_BC = (yc - yb) / (xc - xb);
    double b_BC = yb - a_BC * xb;
    double a_AC = (yc - ya) / (xc - xa);
    double b_AC = ya - a_AC * xa;
    
    printf("Coeficientes:\n");
    printf("AB: a=%.6f, b=%.6f\n", a_AB, b_AB);
    printf("BC: a=%.6f, b=%.6f\n", a_BC, b_BC);
    printf("AC: a=%.6f, b=%.6f\n", a_AC, b_AC);
    printf("Ponto inicial: (%.6f, %.6f)\n", x, y);
    printf("Ângulo inicial: %.1f graus\n", angulo_inicial);
    
    // Número de intersecções
    int max_intersections = 20;
    int intersection_count = 0;
    
    // Nome dos arquivos baseado em N
    char filename[50];
    char filename_x_vx[50];
    sprintf(filename, "CV_N%d.txt", N);
    sprintf(filename_x_vx, "x_vx_N%d.txt", N);
    
    // Arquivos com resultados:
    // CV_N*.txt    -> Guarda a trajetória 2D completa da partícula (X, Y)
    // x_vx_N*.txt  -> Guarda os dados para o mapa de Poincaré (Posição X e Velocidade VX) 
    //                 apenas quando a partícula colide com a base (lado AB).
    FILE *file = fopen(filename, "w");
    FILE *file_x_vx = fopen(filename_x_vx, "w");
    
    if (file == NULL || file_x_vx == NULL) {
        return 1;
    }
    
    // Salvar os pontos dos lados do triângulo primeiro (sem comentários)
    fprintf(file, "%.6f\t%.6f\n", xa, ya);    // Ponto A
    fprintf(file, "%.6f\t%.6f\n", xb, yb);    // Ponto B
    fprintf(file, "%.6f\t%.6f\n", xb, yb);    // Ponto B (para conectar com C)
    fprintf(file, "%.6f\t%.6f\n", xc, yc);    // Ponto C
    fprintf(file, "%.6f\t%.6f\n", xc, yc);    // Ponto C (para conectar com A)
    fprintf(file, "%.6f\t%.6f\n", xa, ya);    // Ponto A (para fechar o triângulo)
    
    // Ponto inicial da trajetória (na borda de baixo)
    fprintf(file, "%.6f\t%.6f\n", x, y);
    
    // Salvar ponto inicial no arquivo x_vx (se estiver no lado AB)
    if (ponto_no_segmento(x, y, xa, ya, xb, yb)) {
        fprintf(file_x_vx, "%.6f\t%.6f\n", x, vx);
    }
    
    // ========================================================================
    // [3] DINÂMICA DO SISTEMA (INTERSEÇÕES E REFLEXÕES)
    // A partir deste ponto ocorre o loop principal que calcula onde a 
    // partícula colide e como ela reflete dentro do triângulo (bilhar).
    // ========================================================================
    while (intersection_count < max_intersections) {
        // Interseção com o lado AB
        double xAB = (b_AB - b) / (a - a_AB);
        double yAB = a * xAB + b;
        
        // Interseção com o lado BC
        double xBC = (b_BC - b) / (a - a_BC);
        double yBC = a * xBC + b;
        
        // Interseção com o lado AC
        double xAC = (b_AC - b) / (a - a_AC);
        double yAC = a * xAC + b;
        
        bool found_intersection = false;
        
        // Verificar interseção com AB (lado de baixo)
        if (ponto_no_segmento(xAB, yAB, xa, ya, xb, yb) && 
            intersecao_valida(xAB, yAB, x, y)) {
            
            printf("Interseção %d com AB: (%.6f, %.6f)\n", intersection_count + 1, xAB, yAB);
            fprintf(file, "%.6f\t%.6f\n", xAB, yAB);
            
            // Salvar no arquivo x_vx apenas para colisões no lado AB
            fprintf(file_x_vx, "%.6f\t%.6f\n", xAB, vx);
            
            x = xAB;
            y = yAB;
            calcular_reflexao(&vx, &vy, a_AB, b_AB);
            a = vy / vx;
            b = y - a * x;
            found_intersection = true;
            
        } 
        // Verificar interseção com BC
        else if (ponto_no_segmento(xBC, yBC, xb, yb, xc, yc) && 
                 intersecao_valida(xBC, yBC, x, y)) {
            
            printf("Interseção %d com BC: (%.6f, %.6f)\n", intersection_count + 1, xBC, yBC);
            fprintf(file, "%.6f\t%.6f\n", xBC, yBC);
            
            x = xBC;
            y = yBC;
            calcular_reflexao(&vx, &vy, a_BC, b_BC);
            a = vy / vx;
            b = y - a * x;
            found_intersection = true;
            
        } 
        // Verificar interseção com AC
        else if (ponto_no_segmento(xAC, yAC, xa, ya, xc, yc) && 
                 intersecao_valida(xAC, yAC, x, y)) {
            
            printf("Interseção %d com AC: (%.6f, %.6f)\n", intersection_count + 1, xAC, yAC);
            fprintf(file, "%.6f\t%.6f\n", xAC, yAC);
            
            x = xAC;
            y = yAC;
            calcular_reflexao(&vx, &vy, a_AC, b_AC);
            a = vy / vx;
            b = y - a * x;
            found_intersection = true;
        }
        
        if (!found_intersection) {
            printf("Erro na iteração %d - posição: (%.6f, %.6f)\n", intersection_count + 1, x, y);
            break;
        }
        
        intersection_count++;
    }
    
    fclose(file);
    fclose(file_x_vx);
    printf("Simulação concluída com %d interseções\n", intersection_count);
    printf("Arquivo da trajetória salvo como: %s\n", filename);
    printf("Arquivo x_vx (apenas lado AB) salvo como: %s\n", filename_x_vx);
    return 0;
}
