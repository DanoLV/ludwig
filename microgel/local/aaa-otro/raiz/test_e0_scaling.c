// Test para entender el escalamiento de e0
// epsilon * e0 * L se agrega como carga ficticia
// El RHS tiene factor eunit*beta
// El potencial resuelto es adimensional

#include <stdio.h>

int main() {
    double e0 = 0.1;           // Campo externo configurado
    double epsilon = 1.0e4;     // Permitividad
    int L = 32;                 // Tamaño del dominio
    double kt = 0.00001;        // Temperatura
    double eunit = 1.0;         // Unidad de carga (default)
    double beta = 1.0 / kt;     // Factor de Boltzmann
    
    // Carga ficticia agregada en el borde
    double rho_ficticia = epsilon * e0 * L;
    
    printf("Configuración:\n");
    printf("  e0 = %.6f\n", e0);
    printf("  epsilon = %.1e\n", epsilon);
    printf("  L = %d\n", L);
    printf("  kT = %.1e\n", kt);
    printf("  beta = %.1e\n", beta);
    printf("  eunit = %.1f\n", eunit);
    printf("\n");
    
    printf("Carga ficticia agregada: rho = epsilon * e0 * L = %.1e\n", rho_ficticia);
    printf("\n");
    
    // El RHS se multiplica por eunit*beta antes del solver
    double rho_adim = rho_ficticia * eunit * beta;
    printf("RHS adimensional: rho * eunit * beta = %.1e\n", rho_adim);
    printf("\n");
    
    // El potencial resultante para campo uniforme sería psi = -e0 * x
    // En unidades adimensionales: psi_adim = psi * eunit * beta / ???
    
    printf("Análisis:\n");
    printf("Si e0 está en unidades de campo 'físico', el campo medido debería ser:\n");
    printf("  E_adim = e0 / kT = %.1e\n", e0 / kt);
    printf("\n");
    printf("Si e0 está en unidades de gradiente adimensional, el campo medido debería ser:\n");
    printf("  E_adim = e0 = %.6f\n", e0);
    printf("\n");
    printf("Campo medido real: ~0.001 - 0.008\n");
    printf("Ratio medido/esperado: %.2f%% - %.2f%%\n", 
           0.001/e0*100, 0.008/e0*100);
    
    return 0;
}
