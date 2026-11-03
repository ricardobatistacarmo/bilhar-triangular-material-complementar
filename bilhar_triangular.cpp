#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

const double PI = 3.14159265358979323846;
const double WIR = 1e-9;

typedef struct {
    double x1, y1;
    double x2, y2;
    const char *nome;
    bool base;
} Aresta;

double produto_vetorial(double ax, double ay, double bx, double by) {
    return ax * by - ay * bx;
}

/*
 * InterseÃ§Ã£o entre o raio P + tV e o segmento Q + uS.
 * HÃ¡ uma colisÃ£o vÃ¡lida quando t > 0 e 0 <= u <= 1.
 */
bool intersecao_raio_segmento(
    double px, double py,
    double vx, double vy,
    const Aresta *aresta,
    double *t, double *u,
    double *xi, double *yi
) {
    double sx = aresta->x2 - aresta->x1;
    double sy = aresta->y2 - aresta->y1;
    double qpx = aresta->x1 - px;
    double qpy = aresta->y1 - py;
    double denominador = produto_vetorial(vx, vy, sx, sy);

    /* Raio e aresta paralelos ou praticamente paralelos. */
    if (fabs(denominador) < WIR) {
        return false;
    }

    *t = produto_vetorial(qpx, qpy, sx, sy) / denominador;
    *u = produto_vetorial(qpx, qpy, vx, vy) / denominador;

    /* t > WIR evita detectar novamente o ponto da colisÃ£o anterior. */
    if (*t <= WIR || *u < -WIR || *u > 1.0 + WIR) {
        return false;
    }

    *xi = px + (*t) * vx;
    *yi = py + (*t) * vy;
    return true;
}

/* ReflexÃ£o especular: v' = v - 2(v.n)n. */
void calcular_reflexao(double *vx, double *vy, const Aresta *aresta) {
    double sx = aresta->x2 - aresta->x1;
    double sy = aresta->y2 - aresta->y1;
    double comprimento = hypot(sx, sy);

    /* Um dos dois vetores normais unitÃ¡rios Ã  aresta. */
    double nx = -sy / comprimento;
    double ny =  sx / comprimento;
    double produto_escalar = (*vx) * nx + (*vy) * ny;

    *vx -= 2.0 * produto_escalar * nx;
    *vy -= 2.0 * produto_escalar * ny;

    /* Corrige pequenas alteraÃ§Ãµes numÃ©ricas no mÃ³dulo da velocidade. */
    double modulo = hypot(*vx, *vy);
    *vx /= modulo;
    *vy /= modulo;
}

int main(void) {
    /* [1] Geometria do bilhar: AC=N, BC=N+1 e AB=N+2. */
    int N = 7;
    double lado_AC = N;
    double lado_BC = N + 1;
    double lado_AB = N + 2;

    double xa = 0.0, ya = 0.0;
    double xb = lado_AB, yb = 0.0;

    /* Coordenadas de C obtidas diretamente das trÃªs medidas dos lados. */
    double xc = (lado_AB * lado_AB + lado_AC * lado_AC
                - lado_BC * lado_BC) / (2.0 * lado_AB);
    double yc_quadrado = lado_AC * lado_AC - xc * xc;

    if (yc_quadrado <= 0.0) {
        fprintf(stderr, "Erro: os lados nÃ£o formam um triÃ¢ngulo vÃ¡lido.\n");
        return 1;
    }

    double yc = sqrt(yc_quadrado);

    Aresta arestas[3] = {
        {xa, ya, xb, yb, "AB", true},
        {xb, yb, xc, yc, "BC", false},
        {xc, yc, xa, ya, "CA", false}
    };

    printf("TriÃ¢ngulo com lados: %.0f, %.0f, %.0f\n",
           lado_AC, lado_BC, lado_AB);
    printf("VÃ©rtices:\n");
    printf("A: (%.6f, %.6f)\n", xa, ya);
    printf("B: (%.6f, %.6f)\n", xb, yb);
    printf("C: (%.6f, %.6f)\n", xc, yc);

    /* [2] CondiÃ§Ãµes iniciais. O Ã¢ngulo Ã© medido a partir do eixo x positivo. */
    double ponto_inicial_x = lado_AB / 2.0;
    double ponto_inicial_y = 0.0;
    double angulo_inicial = 90.0;

    double angulo_rad = angulo_inicial * PI / 180.0;
    double vx = cos(angulo_rad);
    double vy = sin(angulo_rad);
    double x = ponto_inicial_x;
    double y = ponto_inicial_y;

    printf("Ponto inicial: (%.6f, %.6f)\n", x, y);
    printf("Ã‚ngulo inicial: %.6f graus\n", angulo_inicial);

    int max_intersections = 20;
    int intersection_count = 0;

    char filename[50];
    char filename_x_vx[50];
    snprintf(filename, sizeof(filename), "CV_N%d.txt", N);
    snprintf(filename_x_vx, sizeof(filename_x_vx), "x_vx_N%d.txt", N);

    FILE *file = fopen(filename, "w");
    FILE *file_x_vx = fopen(filename_x_vx, "w");

    if (file == NULL || file_x_vx == NULL) {
        fprintf(stderr, "Erro ao abrir os arquivos de saÃ­da.\n");
        if (file != NULL) fclose(file);
        if (file_x_vx != NULL) fclose(file_x_vx);
        return 1;
    }

    /* Segmentos que desenham o contorno do triÃ¢ngulo. */
    fprintf(file, "%.6f\t%.6f\n", xa, ya);
    fprintf(file, "%.6f\t%.6f\n", xb, yb);
    fprintf(file, "%.6f\t%.6f\n", xb, yb);
    fprintf(file, "%.6f\t%.6f\n", xc, yc);
    fprintf(file, "%.6f\t%.6f\n", xc, yc);
    fprintf(file, "%.6f\t%.6f\n", xa, ya);

    fprintf(file, "%.6f\t%.6f\n", x, y);

    /* Registra a condiÃ§Ã£o inicial se ela estiver sobre a base. */
    if (fabs(y - ya) <= WIR && x >= xa - WIR && x <= xb + WIR) {
        fprintf(file_x_vx, "%.6f\t%.6f\n", x, vx);
    }

    /* [3] DinÃ¢mica: encontra sempre a colisÃ£o futura mais prÃ³xima. */
    while (intersection_count < max_intersections) {
        int indice_colisao = -1;
        double menor_t = DBL_MAX;
        double x_colisao = 0.0;
        double y_colisao = 0.0;

        for (int i = 0; i < 3; i++) {
            double t, u, xi, yi;

            if (intersecao_raio_segmento(
                    x, y, vx, vy, &arestas[i],
                    &t, &u, &xi, &yi) && t < menor_t) {
                menor_t = t;
                indice_colisao = i;
                x_colisao = xi;
                y_colisao = yi;
            }
        }

        if (indice_colisao < 0) {
            fprintf(stderr,
                    "Nenhuma colisÃ£o futura encontrada na iteraÃ§Ã£o %d. "
                    "Verifique se o ponto estÃ¡ no interior ou na borda e se "
                    "a velocidade aponta para dentro do triÃ¢ngulo.\n",
                    intersection_count + 1);
            break;
        }

        Aresta *aresta = &arestas[indice_colisao];

        printf("InterseÃ§Ã£o %d com %s: (%.6f, %.6f)\n",
               intersection_count + 1, aresta->nome,
               x_colisao, y_colisao);
        fprintf(file, "%.6f\t%.6f\n", x_colisao, y_colisao);

        /* Para o mapa, registra a posiÃ§Ã£o e vx imediatamente antes da reflexÃ£o. */
        if (aresta->base) {
            fprintf(file_x_vx, "%.6f\t%.6f\n", x_colisao, vx);
        }

        x = x_colisao;
        y = y_colisao;
        calcular_reflexao(&vx, &vy, aresta);
        intersection_count++;
    }

    fclose(file);
    fclose(file_x_vx);

    printf("SimulaÃ§Ã£o concluÃ­da com %d interseÃ§Ãµes.\n", intersection_count);
    printf("TrajetÃ³ria salva em: %s\n", filename);
    printf("SeÃ§Ã£o de PoincarÃ© salva em: %s\n", filename_x_vx);
    return 0;
}
