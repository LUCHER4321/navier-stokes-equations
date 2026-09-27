# Navier Stokes Equations

This repository contains a purely mathematical, highly generalized C# implementation of an incompressible fluid simulation using the Navier-Stokes equations.

The core feature of this simulator is its **$n$-dimensional capability**. Unlike standard grid-based fluid solvers that are hardcoded for 2D or 3D spaces, this solver dynamically adapts to any number of spatial dimensions based entirely on the shape of the initial configuration parameters.

## Mathematical Foundations

The simulation calculates the fluid's acceleration by evaluating the fundamental Navier-Stokes momentum equation. The gradients are computed numerically using the central finite difference method across the $n$-dimensional grid.

### Finite Difference Approximations

The spatial derivatives for pressure gradients, velocity gradients, and the second derivatives of velocity are approximated as follows:

**Pressure Gradient:**

$$
p_{;i}\left( \overrightarrow{x} \right) \approx \frac{p\left( \overrightarrow{x} + \Delta x\overrightarrow{e_i} \right) - p\left( \overrightarrow{x} - \Delta x\overrightarrow{e_i} \right)}{2\Delta x}
$$

**Velocity Gradient (First Derivative):**

$$
\dot{x^i}_{;j}\left( \overrightarrow{x} \right) \approx \frac{\dot{x^i}\left( \overrightarrow{x} + \Delta x\overrightarrow{e_j} \right) - \dot{x^i}\left( \overrightarrow{x} - \Delta x\overrightarrow{e_j} \right)}{2\Delta x}
$$

**Velocity Gradient (Second Derivative):**

$$
\dot{x^i}_{;jk}\left( \overrightarrow{x} \right) \approx \frac{\dot{x^i}\left( \overrightarrow{x} + \Delta x\overrightarrow{e_j} + \Delta x\overrightarrow{e_k} \right) - \dot{x^i}\left( \overrightarrow{x} - \Delta x\overrightarrow{e_j} + \Delta x\overrightarrow{e_k} \right) - \dot{x^i}\left( \overrightarrow{x} + \Delta x\overrightarrow{e_j} - \Delta x\overrightarrow{e_k} \right) + \dot{x^i}\left( \overrightarrow{x} - \Delta x\overrightarrow{e_j} - \Delta x\overrightarrow{e_k} \right)}{4\left( \Delta x \right)^2}
$$

### The Navier-Stokes Momentum Equation

The standard momentum equation for incompressible flow is evaluated tensorially:

$$
\rho\left( \ddot{x^i} + \dot{x^j}\dot{x^i}_{;j} \right) = -g^{ij}p_{;j} + \mu g^{jk}\dot{x^i}_{;jk} + f^i
$$

To integrate the system over time, we isolate the acceleration ($\ddot{x^i}$):

$$
\ddot{x^i} = \frac{-g^{ij}p_{;j} + \mu g^{jk}\dot{x^i}_{;jk} + f^i}{\rho} - \dot{x^j}\dot{x^i}_{;j}
$$

_Where:_

- $\rho$ = Density
- $\mu$ = Dynamic viscosity
- $p$ = Pressure field
- $f^i$ = External forces
- $g^{ij}$ = Metric tensor (Identity matrix in this standard Cartesian space)

## N-Dimensional Architecture

To overcome the static-rank limitations of native multi-dimensional arrays in C# (e.g., `double[,,]`), the simulator utilizes a **Flat Grid Architecture**.

1. **Memory Flattening:** The entire $n$-dimensional grid is flattened into single-dimensional arrays (e.g., `double[] pressureField`, `bool[] fluid`)
2. **Dynamic Tensor Resolution:** Multidimensional arrays combined with jagged arrays (`double[][,]`, `double[][,,]`) as tensors that adapt precisely to the number of dimensions provided during instantiation
3. **Coordinate Mapping:** A generic indexing algorithm automatically maps mathematical coordinates ($x, y, z, w...$) to the 1D flat array memory space

## Usage / Quick Start

The class `IncompressibleFluid` configures its dimensions based entirely on the length of the `extension` array passed to the constructor.

### Example: Setting up a 2D Fluid Simulation

```csharp
using System;
using NavierStokes;

class Program
{
    static void Main()
    {
        // 1. Define dimensions by array length (2 items = 2D)
        double[] dimensions = new double[] { 10.0, 10.0 }; // 10x10 space
        double deltaX = 1.0;     // Grid resolution
        double density = 1.225;  // Fluid density (e.g., air)
        double viscosity = 1.8e-5;
        // 2. Define domain boundary (rectangular/square in this case)
        Func<double[], bool> domainShape = (pos) => true;
        // 3. Initialize Simulator
        IncompressibleFluid sim = new IncompressibleFluid(
            dimensions, deltaX, density, viscosity, domainShape
        );
        // 4. Inject forces or pressure
        double[] forceDir = { 5.0, 0.0 }; // Force along the X-axis
        double[][] positions = { new double[] { 5.0, 5.0 } }; // Applied at the center
        sim.AddForce(forceDir, positions);
        // 5. Advance simulation over time
        double deltaTime = 0.016; // 60 FPS
        sim.Step(deltaTime);
        Console.WriteLine("Simulation step completed.");
    }
}
```

To create a **3D**, **4D**, or higher-dimensional simulation, simply append more elements to the `dimensions` array (e.g., `new double[] { 10.0, 10.0, 10.0 }` for 3D), and ensure your coordinate inputs match the target dimensionality.

## API Reference

### `IncompressibleFluid(double[] extension, double deltaX, double density, double viscosity, Func<double[], bool> shape)`

Constructs the simulation grid. The length of `extension` defines the simulation's dimensionality ($n$). The `shape` delegate evaluates physical coordinates to determine if a grid cell contains fluid or acts as a solid boundary.

### `void Step(double deltaTime)`

Computes the gradients, evaluates the momentum equation, and updates the velocity field across the simulation using an explicit Euler integration step over `deltaTime`.

### `void AddForce(double[] force, double[][] positions)`

Injects external continuous forces at specific physical coordinates.

### `void AddPressure(double pressure, double[][] positions)`

Locally modifies the pressure field at the specified physical coordinates.
