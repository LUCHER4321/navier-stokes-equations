using System;

namespace NavierStokes
{
    public class IncompressibleFluid
    {
        public int dimensions { get; private set; }
        public double deltaX { get; private set; }
        public double[][] velocityField { get; private set; }
        public double[] pressureField;
        public double[][] externalForces;
        public bool[] fluid { get; private set; }
        public double density { get; private set; }
        public double viscosity { get; private set; }
        public int[] sizes { get; private set; }
        public int totalCells { get; private set; }

        public IncompressibleFluid(double[] extension, double deltaX, double density, double viscosity, Func<double[], bool> shape, double[][] initialVelocityField = null)
        {
            this.dimensions = extension.Length;
            this.deltaX = deltaX;
            this.density = density;
            this.viscosity = viscosity;
            this.sizes = new int[dimensions];
            this.totalCells = 1;
            for (int i = 0; i < dimensions; i++)
            {
                sizes[i] = (int)(extension[i] / deltaX);
                totalCells *= sizes[i];
            }
            velocityField = new double[totalCells][];
            pressureField = new double[totalCells];
            externalForces = new double[totalCells][];
            fluid = new bool[totalCells];
            for (int i = 0; i < totalCells; i++)
            {
                int[] coords = GetCoords(i);
                double[] physicalPos = new double[dimensions];
                for (int d = 0; d < dimensions; d++) physicalPos[d] = coords[d] * deltaX;
                velocityField[i] = initialVelocityField == null ? new double[dimensions] : initialVelocityField[i] ?? new double[dimensions];
                externalForces[i] = new double[dimensions];
                fluid[i] = shape(physicalPos);
            }
        }

        public int[] GetCoords(int index)
        {
            int[] coords = new int[dimensions];
            int current = index;
            for (int i = dimensions - 1; i >= 0; i--)
            {
                coords[i] = current % sizes[i];
                current /= sizes[i];
            }
            return coords;
        }

        public int GetIndex(int[] coords)
        {
            int index = 0;
            int multiplier = 1;
            for (int i = dimensions - 1; i >= 0; i--)
            {
                if (coords[i] < 0 || coords[i] >= sizes[i]) return -1;
                index += coords[i] * multiplier;
                multiplier *= sizes[i];
            }
            return index;
        }

        private int GetNeighborIndex(int[] coords, int dim, int offset)
        {
            int original = coords[dim];
            coords[dim] += offset;
            int idx = GetIndex(coords);
            coords[dim] = original;
            return idx;
        }

        private double Derivative(Func<int, double> getValue, int centerIdx, int[] coords, int dim)
        {
            int nextIdx = GetNeighborIndex(coords, dim, 1);
            int prevIdx = GetNeighborIndex(coords, dim, -1);
            bool hasNext = nextIdx != -1 && fluid[nextIdx];
            bool hasPrev = prevIdx != -1 && fluid[prevIdx];
            if (hasNext && hasPrev) return (getValue(nextIdx) - getValue(prevIdx)) / (2.0 * deltaX);
            else if (hasNext) return (getValue(nextIdx) - getValue(centerIdx)) / deltaX;
            else if (hasPrev) return (getValue(centerIdx) - getValue(prevIdx)) / deltaX;
            return 0.0;
        }

        private double[][] PressureGradient()
        {
            double[][] gradient = new double[totalCells][];
            for (int i = 0; i < totalCells; i++)
            {
                gradient[i] = new double[dimensions];
                if (!fluid[i]) continue;
                int[] coords = GetCoords(i);
                for (int d = 0; d < dimensions; d++) gradient[i][d] = Derivative(idx => pressureField[idx], i, coords, d);
            }
            return gradient;
        }

        private double[][,] VelocityGradient()
        {
            double[][,] gradient = new double[totalCells][,];
            for (int i = 0; i < totalCells; i++)
            {
                gradient[i] = new double[dimensions, dimensions];
                if (!fluid[i]) continue;
                int[] coords = GetCoords(i);
                for (int spatialDim = 0; spatialDim < dimensions; spatialDim++) for (int velDim = 0; velDim < dimensions; velDim++) gradient[i][spatialDim, velDim] = Derivative(idx => velocityField[idx][velDim], i, coords, spatialDim);
            }
            return gradient;
        }

        private double[][,,] VelocitySecondGradient()
        {
            double[][,,] gradient = new double[totalCells][,,];
            double[][,] firstGradient = VelocityGradient();
            for (int i = 0; i < totalCells; i++)
            {
                gradient[i] = new double[dimensions, dimensions, dimensions];
                if (!fluid[i]) continue;
                int[] coords = GetCoords(i);
                gradient[i] = new double[dimensions, dimensions, dimensions];
                for (int spatialDim1 = 0; spatialDim1 < dimensions; spatialDim1++) for (int velDim = 0; velDim < dimensions; velDim++) for (int spatialDim2 = 0; spatialDim2 < dimensions; spatialDim2++) gradient[i][spatialDim1, velDim, spatialDim2] = Derivative(idx => firstGradient[idx][spatialDim1, velDim], i, coords, spatialDim2);
            }
            return gradient;
        }

        private double[][] Acceleration()
        {
            double[][] acc = new double[totalCells][];
            double[][,] velGrad = VelocityGradient();
            double[][] presGrad = PressureGradient();
            double[][,,] velSecondGrad = VelocitySecondGradient();
            for (int i = 0; i < totalCells; i++)
            {
                acc[i] = new double[dimensions];
                if (!fluid[i]) continue;
                double[] velLaplacian = new double[dimensions];
                for (int velDim = 0; velDim < dimensions; velDim++)
                {
                    double sum = 0;
                    for (int spatialDim = 0; spatialDim < dimensions; spatialDim++) sum += velSecondGrad[i][spatialDim, velDim, spatialDim];
                    velLaplacian[velDim] = sum;
                }
                double[] velDotVelGrad = new double[dimensions];
                for (int velDim = 0; velDim < dimensions; velDim++)
                {
                    double sum = 0;
                    for (int spatialDim = 0; spatialDim < dimensions; spatialDim++) sum += velocityField[i][spatialDim] * velGrad[i][spatialDim, velDim];
                    velDotVelGrad[velDim] = sum;
                }
                for (int d = 0; d < dimensions; d++) acc[i][d] = (-presGrad[i][d] + viscosity * velLaplacian[d] + externalForces[i][d]) / density - velDotVelGrad[d];
            }
            return acc;
        }

        public void Step(double deltaTime)
        {
            double[][] acc = Acceleration();
            for (int i = 0; i < totalCells; i++)
            {
                if (!fluid[i]) continue;
                for (int d = 0; d < dimensions; d++) velocityField[i][d] += acc[i][d] * deltaTime;
            }
        }

        public void AddPressure(double pressure, double[][] positions)
        {
            foreach (double[] pos in positions)
            {
                int[] coords = new int[dimensions];
                for (int d = 0; d < dimensions; d++) coords[d] = (int)(pos[d] / deltaX);
                int idx = GetIndex(coords);
                if (idx != -1 && fluid[idx]) pressureField[idx] += pressure;
            }
        }

        public void AddForce(double[] force, double[][] positions)
        {
            double forceVolume = 0;
            double cellVolume = Math.Pow(deltaX, dimensions);
            foreach (double[] pos in positions)
            {
                int[] coords = new int[dimensions];
                for (int d = 0; d < dimensions; d++) coords[d] = (int)(pos[d] / deltaX);
                int idx = GetIndex(coords);
                if (idx != -1 && fluid[idx]) forceVolume += cellVolume;
            }
            foreach (double[] pos in positions)
            {
                int[] coords = new int[dimensions];
                for (int d = 0; d < dimensions; d++) coords[d] = (int)(pos[d] / deltaX);
                int idx = GetIndex(coords);
                if (idx != -1 && fluid[idx]) for (int d = 0; d < dimensions; d++) externalForces[idx][d] += force[d] / forceVolume;
            }
        }
    }
}