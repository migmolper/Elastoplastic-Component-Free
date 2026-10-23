# Elmt15 — Gradient-Enhanced Nonlocal Damage (FEAP)

Elemento de usuario 3D de **deformación finita** para [FEAP](https://projects.ce.berkeley.edu/feap/) que implementa un modelo de **daño no local con mejora por gradiente** (*Gradient-Enhanced Damage*, GED).

Archivo principal: [`elmt15.f`](elmt15.f)

---

## Descripción

El elemento resuelve un sólido hiperelástico dañado acoplado a un campo no local \(\phi\). Combina:

- Elasticidad **Neo-Hookean** compressible
- Función de daño exponencial \(f_d(K)\)
- Regularización no local mediante \(\phi\) y \(\nabla\phi\) (modelo GED)

Soporta dos formulaciones cinemáticas seleccionables por el parámetro `etype`:

| `etype` | Subrutina      | Formulación              | Tensiones / módulos        |
|---------|----------------|--------------------------|----------------------------|
| `1`     | `GEDmodel`     | Updated Lagrangian (UL)  | Cauchy / espacial         |
| `2`     | `GEDmodelTL`   | Total Lagrangian (TL)    | 2º Piola–Kirchhoff / material |

> **Nota:** comentarios antiguos en el código etiquetan `etype` como “GED vs estándar”; en la implementación actual **ambos caminos son GED**.

---

## Requisitos de malla / control FEAP

| Parámetro | Valor recomendado | Comentario |
|-----------|-------------------|------------|
| `ndm`     | `3`               | Problema 3D |
| `ndf`     | `4`               | 3 desplazamientos + 1 \(\phi\) no local |
| `nen`     | `8`               | Hexaedro de 8 nodos |
| Cuadratura | \(2\times2\times2\) | 8 puntos de Gauss |

Grados de libertad por nodo:

\[
\mathbf{u}_I = (u_x,\, u_y,\, u_z,\, \phi)_I
\]

El ensamblado intercala DOFs como `[ux, uy, uz, φ]` en cada nodo.

---

## Parámetros de material (`isw = 1`)

En el bloque `mate` de FEAP:

```text
mate 1
  user 15
  E, nu, Kd, ndp, etype
```

| Campo   | Significado |
|---------|-------------|
| `E`     | Módulo de Young |
| `nu`    | Coeficiente de Poisson |
| `Kd`    | Umbral de daño (almacenado en `d(243)`) |
| `ndp`   | Parámetro de la ley exponencial de daño (`d(244)`) |
| `etype` | `1` = UL, `2` = TL (`d(245)`) |

A partir de \(E\) y \(\nu\) se obtienen:

\[
\mu = \frac{E}{2(1+\nu)},\qquad
\lambda = \frac{E}{3(1-2\nu)} - \frac{2}{3}\mu
\]

almacenados en `d(22)` (\(\mu\)) y `d(21)` (\(\lambda\)).

### Ejemplo (`Iex.txt`)

```text
0,0,0,3,4,8

parameter
  e  = 440.5d0
  nu = 0.25
  np = 1.
  kd = 1.0

mate 1
  user 15
  e,nu,kd,np,1
```

Conviene inicializar el DOF de \(\phi\) a \(K_d\) (véase `INITial DISP` en `Iex.txt`).

---

## Modelo constitutivo

### Neo-Hookean

Energía elástica (esquema UL, con \(\mathbf{b}\)):

\[
\psi_E = \frac{\mu}{2}(\mathrm{tr}\,\mathbf{b}-3)
        - \mu\ln J
        + \frac{\lambda}{2}(\ln J)^2
\]

- **UL** (`StandardNeoh3`): tensión de Cauchy  
  \(\boldsymbol{\sigma}_E = \mu\,\mathbf{b}/J + (\lambda\ln J - \mu)\,J^{-1}\mathbf{I}\)
- **TL** (`StandardNeoh3reference`): 2º Piola–Kirchhoff  
  \(\mathbf{S}_E = \mu\,\mathbf{I} + (\lambda\ln J - \mu)\,\mathbf{C}^{-1}\)

### Función de daño

```text
fd(K)  = exp(ndp * (Kd - K))
fd'(K) = -ndp * fd(K)
```

Medida de daño proyectada: \(D = 1 - f_d\).

Tensión dañada (con regularización numérica):

\[
\boldsymbol{\sigma}_D = (f_d + 10^{-4})\,\boldsymbol{\sigma}_E
\]

### Evolución local de \(K\)

Por punto de Gauss se almacena la variable interna \(K\) (historia `nhv = 1`), inicializada a \(K_d\).

Parámetros no locales (hardcoded en `LocalNewtonDamage*`):

| Símbolo   | Valor típico en código | Rol |
|-----------|------------------------|-----|
| \(c_d\)   | `10`                   | Longitud / peso del gradiente |
| \(\gamma_d\) | `1`                 | Acoplamiento local–no local |
| \(\beta_d\)  | `1000`              | Penalización |

Fuerza conjugada y superficie de daño (esquema UL):

\[
q = \psi_E + \frac{\beta_d\gamma_d}{n_{dp} f_d}(\phi - \gamma_d K),\qquad
\chi = q - K
\]

Si \(\chi > 0\) (y \(f_d\) no es demasiado pequeña), se actualiza \(K_{n+1}\) con un **Newton local**.

---

## Estructura del código

| Subrutina | Rol |
|-----------|-----|
| `elmt15` | Driver FEAP (`isw`); lee propiedades y despacha UL/TL |
| `GEDmodel` | Residuo y tangente — Updated Lagrangian |
| `GEDmodelTL` | Residuo y tangente — Total Lagrangian |
| `StandardNeoh3` / `StandardNeoh3reference` | Neo-Hookean (Cauchy / 2º PK) |
| `LocalNewtonDamage` / `LocalNewtonDamageReference` | Actualización de daño + derivadas |
| `damageFunction` | \(f_d(K)\) y \(f_d'\) |
| `getKuphi`, `getKphiphi`, `getFphi*` | Acoplamientos \(K_{u\phi}\), \(K_{\phi\phi}\) y residuo en \(\phi\) |
| `getKuu_Fu`, `getKuuGeom` | Parte material y geométrica de \(K_{uu}\) |
| `getShape3D8`, `getDerivativeShape3D8`, `quadrature3D8` | Forma, derivadas y cuadratura del hex8 |
| `getDeformGrad3D`, `getBmat3D`, `getBmat3DNL` | Cinemática y matrices B |
| `gaussValueProjection` | Proyección a nodos (`isw = 8`) |

Variables de historia: `nh1 = 8` (un \(K\) por Gauss point).  
Campos proyectados (`istv = 8`): componentes de tensión + daño \(D\).

---

## Switches FEAP relevantes (`isw`)

| `isw` | Acción en este elemento |
|-------|-------------------------|
| `0`   | Mensaje de tipo de elemento |
| `1`   | Lectura de materiales, checks `ndf`/`nen`/`etype` |
| `3`, `6` | Residual y tangente (vía `GEDmodel` / `GEDmodelTL`) |
| `4`   | Salida de esfuerzos |
| `8`   | Proyección / plot |
| `14`  | Inicialización de historia |

---

## Archivos del repositorio

| Archivo | Contenido |
|---------|-----------|
| `elmt15.f` | Elemento GED (daño no local + gradiente) |
| `elmt16.f` | Versión simplificada: solo Neo-Hookean (sin daño) |
| `Iex.txt` | Ejemplo de input FEAP (para elmt15) |
| `hc.nodes` / `hc.elems` | Coordenadas y conectividades del ejemplo |

### Elemento simplificado `elmt16` (solo Neo-Hookean)

Copia reducida de `elmt15` con formulación **Updated Lagrangian**, hex8 y Gauss \(2\times2\times2\). Sin daño, sin \(\phi\) no local y sin historia.

```text
0,0,0,3,3,8          ! ndm=3, ndf=3, nen=8

mate 1
  user 16
  E, nu
```

Las rutinas auxiliares llevan prefijo `nh` para poder enlazar `elmt15.f` y `elmt16.f` a la vez.

---

## Historial (del propio código)

| Fecha        | Cambio |
|--------------|--------|
| 01/11/2006   | Versión original |
| 09/09/2014   | Implementación tipo elmt14 / basada en `disp30.f` |
| 13/03/2017   | Modelo GED implementado |

---

## Compilación / uso

1. Incluir `elmt15.f` en la librería de elementos de usuario de FEAP (o enlazarlo según el `makefile` del proyecto).
2. Definir el control con `ndm=3`, `ndf=4`, `nen=8`.
3. Declarar el material con `user 15` y los cinco parámetros indicados arriba.
4. Condicionar el 4º DOF (\(\phi\)) en contornos según el problema (p. ej. `cboun` / `eboun` con cuatro flags).

---

## Avisos

- El check interno solo exige `ndf ≥ 3`, pero el campo \(\phi\) **requiere `ndf ≥ 4`**.
- Helpers EAS (`kEASMlocalGE3d`, `kEASTransfGE3D`) aparecen en el camino TL pero **no se ensamblan** en \(K\)/\(r\).
- Parámetros \(c_d\), \(\gamma_d\), \(\beta_d\) están fijados en el código; no se leen del input FEAP.
