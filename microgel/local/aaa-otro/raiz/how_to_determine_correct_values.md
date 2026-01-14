# ¿Cómo Determinar el Valor Correcto del Laplaciano?

## La Pregunta Fundamental

**¿Cómo sabemos que el Laplaciano debe dar -6.0 y no otro valor?**

## Método: Campo de Prueba con Solución Analítica Conocida

### Paso 1: Elegir un Campo de Prueba Simple

Usamos: **ψ(x,y,z) = x² + y² + z²**

¿Por qué este campo?
- ✓ La solución analítica es trivial de calcular
- ✓ Es simétrico en las tres direcciones
- ✓ Es cuadrático (perfecto para diferencias finitas de 2do orden)

### Paso 2: Calcular la Solución Analítica (Cálculo Diferencial Básico)

#### A. GRADIENTE:
```
∂ψ/∂x = ∂(x² + y² + z²)/∂x = 2x
∂ψ/∂y = ∂(x² + y² + z²)/∂y = 2y
∂ψ/∂z = ∂(x² + y² + z²)/∂z = 2z

Por tanto: ∇ψ = (2x, 2y, 2z)
```

**Verificación en el punto (1,1,1):**
```
∇ψ|(1,1,1) = (2·1, 2·1, 2·1) = (2, 2, 2)  ✓
```

#### B. LAPLACIANO:
```
∂²ψ/∂x² = ∂(2x)/∂x = 2
∂²ψ/∂y² = ∂(2y)/∂y = 2
∂²ψ/∂z² = ∂(2z)/∂z = 2

Por tanto: ∇²ψ = ∂²ψ/∂x² + ∂²ψ/∂y² + ∂²ψ/∂z²
               = 2 + 2 + 2
               = 6
```

**Resultado analítico: ∇²ψ = 6.0** (constante en todo el espacio)

### Paso 3: Verificar con Diferencias Finitas Estándar

Las diferencias finitas centradas de 2do orden con espaciado h=1:

```
∂²ψ/∂x² ≈ [ψ(x+1,y,z) - 2ψ(x,y,z) + ψ(x-1,y,z)] / h²
```

Con h=1:
```
∂²ψ/∂x² ≈ ψ(x+1,y,z) - 2ψ(x,y,z) + ψ(x-1,y,z)
```

**Verificación manual en el origen (0,0,0):**

```
ψ(0,0,0) = 0² + 0² + 0² = 0

Dirección x:
  ψ(1,0,0) = 1² + 0² + 0² = 1
  ψ(-1,0,0) = (-1)² + 0² + 0² = 1
  ∂²ψ/∂x² = 1 - 2(0) + 1 = 2  ✓

Dirección y:
  ψ(0,1,0) = 0² + 1² + 0² = 1
  ψ(0,-1,0) = 0² + (-1)² + 0² = 1
  ∂²ψ/∂y² = 1 - 2(0) + 1 = 2  ✓

Dirección z:
  ψ(0,0,1) = 0² + 0² + 1² = 1
  ψ(0,0,-1) = 0² + 0² + (-1)² = 1
  ∂²ψ/∂z² = 1 - 2(0) + 1 = 2  ✓

Total:
  ∇²ψ = 2 + 2 + 2 = 6  ✓✓✓
```

**Conclusión: Las diferencias finitas confirman ∇²ψ = 6.0**

### Paso 4: Convención de Signos en Ludwig

En Ludwig, los pesos `wlaplacian` incorporan el signo negativo:
```c
// La convención es que wlaplacian da -∇²ψ
resultado_codigo = Σ wlaplacian[p] · ψ[p] = -∇²ψ
```

Por tanto, para ψ = x² + y² + z²:
```
Resultado esperado = -6.0
```

### Paso 5: Verificar con D3Q7 (Stencil de Referencia)

D3Q7 es el stencil más simple y sabemos que está correcto:

```
Puntos del stencil D3Q7:
  p=0: (0, 0, 0) → ψ = 0
  p=1: (1, 0, 0) → ψ = 1
  p=2: (0, 1, 0) → ψ = 1
  p=3: (0, 0, 1) → ψ = 1
  p=4: (0, 0,-1) → ψ = 1
  p=5: (0,-1, 0) → ψ = 1
  p=6: (-1, 0, 0) → ψ = 1

Pesos LB base:
  w[0] = 2/8 = 0.25
  w[1..6] = 1/8 = 0.125

Con α = -8.0:
  wlaplacian[0] = -(-8.0)·6·(1/8) = 6.0
  wlaplacian[1..6] = -8.0·(1/8) = -1.0

Cálculo:
  Σ wlaplacian[p]·ψ[p] = 6.0·0 + 6·(-1.0)·1
                        = -6.0  ✓
```

**D3Q7 con α = -8.0 da -6.0 ✓**

### Paso 6: Aplicar el Mismo Criterio a D3Q19 y D3Q27

**Principio de Consistencia:**
> Todos los stencils deben dar el MISMO resultado para el MISMO campo físico

Si D3Q7 da -6.0 (correcto), entonces D3Q19 y D3Q27 también DEBEN dar -6.0.

#### Para D3Q19:

```
Pesos LB base: w[0]=12/36, w[1]=1/36, w[3]=2/36, etc.

Valores de ψ en el stencil:
  ψ[0] = 0
  ψ en vecinos axiales (ej: [1,0,0]) = 1
  ψ en vecinos diagonales (ej: [1,1,0]) = 2

Σ(p>0) w[p]·ψ[p] = 1.0

Para obtener -6.0:
  α · 1.0 = -6.0
  α = -6.0  ✓
```

**Verificación numérica:**
```python
import numpy as np

cv = [[0,0,0], [1,1,0], [1,0,1], [1,0,0], ...]  # D3Q19
wv = [12/36, 1/36, 1/36, 2/36, ...]

psi = [cv[p][0]**2 + cv[p][1]**2 + cv[p][2]**2 for p in range(19)]
# psi = [0, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2]

wlap = -6.0 * wv
wlap[0] = -sum(wlap[1:])

result = sum(wlap[p] * psi[p] for p in range(19))
# result = -6.0  ✓
```

#### Para D3Q27:

```
Mismo procedimiento:
  Σ(p>0) w[p]·ψ[p] = 1.0

  α · 1.0 = -6.0
  α = -6.0  ✓
```

## Resumen: ¿Cómo se Determina el Valor Correcto?

### Es un Proceso en 6 Pasos:

1. **Elegir campo de prueba**: ψ = x² + y² + z² (solución conocida)

2. **Calcular solución analítica**: ∇²ψ = 6.0 (cálculo diferencial)

3. **Verificar con diferencias finitas**: Confirmar que FD da 6.0

4. **Aplicar convención de signos**: Resultado esperado = -6.0

5. **Verificar stencil de referencia**: D3Q7 con α=-8.0 da -6.0 ✓

6. **Imponer consistencia**: Todos los stencils deben dar -6.0

### NO es una Elección Arbitraria

El valor -6.0 proviene de:
- ✅ **Matemática pura**: Cálculo diferencial de ∂²ψ/∂x²
- ✅ **Diferencias finitas estándar**: Fórmulas bien conocidas
- ✅ **Consistencia física**: Misma ecuación → mismo resultado
- ✅ **Verificación numérica**: Podemos calcularlo explícitamente

### Valores Correctos Derivados:

| Stencil | Factor α (Laplaciano) | Factor β (Gradiente) | Resultado |
|---------|----------------------|---------------------|-----------|
| D3Q7    | -8.0                 | +4.0                | ∇²ψ = -6.0 ✓ |
| D3Q19   | **-6.0** (era -36.0) | +3.0                | ∇²ψ = -6.0 ✓ |
| D3Q27   | **-6.0** (era -216.0)| +3.0                | ∇²ψ = -6.0 ✓ |

## Verificación Adicional: Test de Gradiente

Para completar, verificamos que el gradiente también da el valor correcto:

**Campo de prueba en (1,1,1):**
```
ψ = x² + y² + z²
∇ψ = (2x, 2y, 2z) = (2, 2, 2)
E = -∇ψ = (-2, -2, -2)
```

**Con D3Q19 y β = 3.0:**
```python
wgrad = 3.0 * wv
wgrad[0] = 0.0

# Evaluar en (1,1,1)
psi_111 = [(1+cv[p][0])**2 + (1+cv[p][1])**2 + (1+cv[p][2])**2
           for p in range(19)]

E_x = -sum(wgrad[p] * cv[p][0] * psi_111[p] for p in range(19))
# E_x = -2.0  ✓
```

**Conclusión: El gradiente con β=3.0 ya estaba correcto.**

## Mensaje Final

El valor correcto no es una "opinión" ni una "convención arbitraria". Es el resultado de:

1. Aplicar correctamente la definición matemática del Laplaciano
2. Usar diferencias finitas de 2do orden estándar
3. Imponer consistencia entre diferentes stencils
4. Verificar con cálculos explícitos

**Los valores anteriores (-36.0 y -216.0) estaban matemáticamente incorrectos.**

El fix aplicado corrige este error fundamental y garantiza que todos los stencils den resultados físicamente consistentes.
