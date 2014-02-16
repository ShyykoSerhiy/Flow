package solver;

import draw.Point;
import shape.Shape;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static java.lang.Math.*;

/**
 * Created with IntelliJ IDEA.
 * User: shyyko
 * Date: 08.02.14
 * Time: 15:17
 */
public class Solver {
    private Shape shape;
    /**
     * Is vector which has length 1 i.e. (x^2 +y^2)^(1/2) == 1
     */
    private Point vEternity;

    /**
     * Hamma_j
     */
    private List<Double> hamma;
    private double hamma0;
    private double alpha;

    public Solver(Shape shape) {
        this.shape = shape;
        alpha = 0; //todo
        this.vEternity = new Point(cos(alpha), sin(alpha)); //todo move to constructor parameters

        hamma = new ArrayList<Double>(Collections.nCopies(shape.getKolokPoints().size(), 0.0));
        hamma0 = 10; //todo move to constructor parameters
    }

    public void solve() {
        List<List<Double>> matrix = new ArrayList<List<Double>>();
        List<Double> rightPartOfSystem = new ArrayList<Double>();
        int amountOfCollocPoints = shape.getKolokPoints().size();

        for (int i = 0; i < amountOfCollocPoints; i++) {
            Point collocPoint = shape.getKolokPoints().get(i);
            List<Double> rowOfMatrix = new ArrayList<Double>(amountOfCollocPoints);
            for (int j = 0; j < amountOfCollocPoints + 1; j++) {
                Point vJ = computeVj(collocPoint, shape.getAllPoints().get(j));
                double coefficientToAdd = vJ.multiply(shape.getNormal().get(i));
                rowOfMatrix.add(coefficientToAdd);
            }
            matrix.add(rowOfMatrix);
            rightPartOfSystem.add(-vEternity.multiply(shape.getNormal().get(i)));
        }
        matrix.add(new ArrayList<Double>(Collections.nCopies(amountOfCollocPoints + 1, 1.0)));
        rightPartOfSystem.add(hamma0);

        hamma = gaussMethod(matrix, rightPartOfSystem);
        System.out.println(hamma);
    }

    public double getPhi(Point point) {
        double toReturn = point.getX() * cos(alpha) + point.getY() * sin(alpha);

        double sum = 0;
        for (int i = 0; i < hamma.size(); i++) {
            sum += hamma.get(i) * atan2(
                    point.getY() - shape.getAllPoints().get(i).getY(),
                    point.getY() - shape.getAllPoints().get(i).getY()
            );
        }
        toReturn += 1 / (2 * PI) * sum;
        return toReturn;
    }

    public double getPsi(Point point) {
        double toReturn = cos(alpha) * point.getY() - sin(alpha) * point.getX();

        double sum = 0;
        for (int i = 0; i < hamma.size(); i++) {
            sum += hamma.get(i) * log(computeRj(point, shape.getAllPoints().get(i)));
        }
        toReturn -= 1 / (2 * PI) * sum;
        return toReturn;
    }

    private Point computeVj(Point firstPoint, Point secondPoint) {
        return new Point(
                getU(firstPoint, secondPoint),
                getV(firstPoint, secondPoint)
        );
    }

    public Point computeV(Point point) {
        Point sum = vEternity.add(new Point(0, 0));
        for (int i = 0; i < hamma.size(); i++) {
            sum = sum.add(computeVj(point, shape.getAllPoints().get(i)).multiply(hamma.get(i)));
        }
        return sum;
    }

    private double getU(Point firstPoint, Point secondPoint) {
        return 1 / (2 * PI) * ((secondPoint.getY() - firstPoint.getY()) / pow(computeRj(firstPoint, secondPoint), 2));
    }

    private double getV(Point firstPoint, Point secondPoint) {
        return 1 / (2 * PI) * ((firstPoint.getX() - secondPoint.getX()) / pow(computeRj(firstPoint, secondPoint), 2));
    }

    private double computeRj(Point firstPoint, Point secondPoint) {
        return max(getSigma(), sqrt(
                pow(firstPoint.getX() - secondPoint.getX(), 2) + pow(firstPoint.getY() - secondPoint.getY(), 2))
        );
    }

    public double getSigma() { //todo implement
        return 0.001;
    }


    private static List<Double> gaussElimination(List<List<Double>> matrix, List<Double> rightSide) {
        // Gaussian elimination with partial pivoting
        int size = rightSide.size();

        for (int p = 0; p < size; p++) {
            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < size; i++) {
                if (abs(matrix.get(i).get(p)) > abs(matrix.get(max).get(p))) {
                    max = i;
                }
            }
            List<Double> temp = matrix.get(p);
            matrix.set(p, matrix.get(max));
            matrix.set(max, temp);
            double t = rightSide.get(p);
            rightSide.set(p, rightSide.get(max));
            rightSide.set(max, t);

            // singular or nearly singular
            if (abs(matrix.get(p).get(p)) <= 1e-10) {
                throw new RuntimeException("Matrix is singular or nearly singular");
            }

            // pivot within A and b
            for (int i = p + 1; i < size; i++) {
                double alpha = matrix.get(i).get(p) / matrix.get(p).get(p);
                rightSide.set(i, rightSide.get(i) - alpha * rightSide.get(p));
                for (int j = p; j < size; j++) {
                    matrix.get(i).set(j, matrix.get(i).get(j) - alpha * matrix.get(p).get(j));
                }
            }
        }

        // back substitution
        Double[] x = new Double[size];
        for (int i = size - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < size; j++) {
                sum += matrix.get(i).get(j) * x[j];
            }
            x[i] = (rightSide.get(i) - sum) / matrix.get(i).get(i);
        }
        return new ArrayList<Double>(Arrays.asList(x));
    }

    public static final float EPS = 0.001f;


    private static List<Double> gaussMethod(List<List<Double>> matrix, List<Double> rightSide) {
        float[][] matrixA = new float[rightSide.size()][rightSide.size()];
        float[] b = new float[rightSide.size()];

        for (int i = 0; i < rightSide.size(); i++) {
            for (int j = 0; j < rightSide.size(); j++) {
                matrixA[i][j] = matrix.get(i).get(j).floatValue();
            }
            b[i] = rightSide.get(i).floatValue();
        }
        gaussMethod(matrixA, b);
        List<Double> toReturn = new ArrayList<Double>(rightSide.size());
        for (float value : b) {
            toReturn.add(Double.valueOf(value));
        }
        return toReturn;
    }

    public static void gaussMethod(float[][] a, float[] b) {
        for (int i = 0; i < b.length - 1; i++) {
            int l = 0;
            while (Math.abs(a[i][i]) < EPS) {
                l++;
                if (i < b.length - l) {
                    for (int j = 0; j < b.length; j++) {
                        a[i][j] += a[i + l][j];
                    }
                    b[i] += b[i + l];
                } else {
                    for (int j = 0; j < b.length; j++) {
                        a[i][j] += a[i - l][j];
                    }
                    b[i] += b[i - l];
                }
            }

            for (int j = i + 1; j < b.length; j++) {
                if (Math.abs(a[j][i]) > EPS) {
                    final float koef = -a[j][i] / a[i][i];
                    for (int k = 0; k < b.length; k++) {
                        a[j][k] += koef * a[i][k];
                    }
                    b[j] += b[i] * koef;
                }
            }
        }
        for (int j = b.length - 1; j >= 1; j--) {
            for (int i = 0; i < j; i++) {
                if (Math.abs(a[i][j]) > EPS) {
                    final float koef = -a[i][j] / a[j][j];
                    a[i][j] += koef * a[j][j];
                    b[i] += b[j] * koef;
                }
            }
            b[j] /= a[j][j];
            a[j][j] = 1;
        }
        b[0] /= a[0][0];
        a[0][0] = 1;
    }
}
