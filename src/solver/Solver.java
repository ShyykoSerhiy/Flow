package solver;

import draw.Point;
import draw.PointWithHamma;
import shape.Shape;

import java.util.*;

import static java.lang.Math.*;

/**
 * Created with IntelliJ IDEA.
 * User: shyyko
 * Date: 08.02.14
 * Time: 15:17
 */
public class Solver {
    public static final double DELTA = 0.001;

    private Shape shape;
    /**
     * Is vector which has length 1 i.e. (x^2 +y^2)^(1/2) == 1
     */
    private Point vEternity;

    /**
     * Points
     */
    private Map<PointWithHamma, List<PointWithHamma>> flyingPointsWithHama;
    private double alpha;

    public Solver(Shape shape) {
        this.shape = shape;
        alpha = 0; //todo
        this.vEternity = new Point(cos(alpha), sin(alpha)); //todo move to constructor parameters

        flyingPointsWithHama = new HashMap<PointWithHamma, List<PointWithHamma>>();
        for (PointWithHamma pointWithHamma : shape.getPointsOfSeparation()) {
            flyingPointsWithHama.put(pointWithHamma, new ArrayList<PointWithHamma>());
        }
    }

    public Map<PointWithHamma, List<PointWithHamma>> getFlyingPointsWithHama() {
        return flyingPointsWithHama;
    }

    public void solve() {
        List<List<Double>> matrix = new ArrayList<List<Double>>();
        List<Double> rightPartOfSystem = new ArrayList<Double>();
        int amountOfCollocPoints = shape.getKolokPoints().size();

        for (int i = 0; i < amountOfCollocPoints; i++) {
            Point kolokPoint = shape.getKolokPoints().get(i);
            Point normalInKolocPoint = shape.getNormal().get(i);

            List<Double> rowOfMatrix = new ArrayList<Double>(amountOfCollocPoints);
            for (int j = 0; j < amountOfCollocPoints + 1; j++) {
                Point vJ = computeVj(kolokPoint, shape.getAllPoints().get(j));
                double coefficientToAdd = vJ.multiply(normalInKolocPoint);
                rowOfMatrix.add(coefficientToAdd);
            }
            matrix.add(rowOfMatrix);

            double rightPart = -vEternity.multiply(shape.getNormal().get(i));
            for (List<PointWithHamma> list : flyingPointsWithHama.values()) {
                for (PointWithHamma pointWithHamma : list) {
                    rightPart -= pointWithHamma.getHamma() * computeVj(kolokPoint, pointWithHamma).multiply(normalInKolocPoint);
                }
            }

            rightPartOfSystem.add(rightPart);
        }

        matrix.add(new ArrayList<Double>(Collections.nCopies(amountOfCollocPoints + 1, 1.0)));
        double lastRightPart = 0;
        for (List<PointWithHamma> list : flyingPointsWithHama.values()) {
            for (PointWithHamma pointWithHamma : list) {
                lastRightPart -= pointWithHamma.getHamma();
            }
        }


        rightPartOfSystem.add(lastRightPart);

        List<Double> hamma = gaussMethod(matrix, rightPartOfSystem);
        for (int i = 0; i < hamma.size(); i++) {
            shape.getAllPoints().get(i).setHamma(hamma.get(i));
        }

        // adding flying points (pointsWithHamma)
        for (PointWithHamma pointWithHamma : shape.getPointsOfSeparation()) {
            PointWithHamma flyingPoint = new PointWithHamma(pointWithHamma.getX(), pointWithHamma.getY(), pointWithHamma.getHamma());
            flyingPointsWithHama.get(pointWithHamma).add(flyingPoint);
        }
    }

    /**
     * Moving flying points with intersection check
     */
    public void movePoints() {
        if (flyingPointsWithHama.get(shape.getPointsOfSeparation().get(0)).isEmpty()) {
            return;
        }

        double maxSpeed = 0;
        for (List<PointWithHamma> flyingPoints : flyingPointsWithHama.values()) {
            PointWithHamma pointWithHamma = flyingPoints.get(flyingPoints.size() - 1);
            Point speed = computeV(pointWithHamma);
            double speedDouble = Point.distance(speed, new Point(0, 0));
            if (maxSpeed < speedDouble) {
                maxSpeed = speedDouble;
            }
        }
        double time = shape.getDelta() / maxSpeed;

        for (List<PointWithHamma> flyingPoints : flyingPointsWithHama.values()) {
            for (int i = 0; i < flyingPoints.size(); i++) {
                PointWithHamma oldPoint = flyingPoints.get(i);
                Point newPointCoordinate = null;
                newPointCoordinate = computeV(oldPoint).multiply(time).add(oldPoint);
                newPointCoordinate = reflectIsCrossCarcases(oldPoint, newPointCoordinate);
                newPointCoordinate = awayFromSegmentationPoint(newPointCoordinate);

                flyingPoints.set(i, new PointWithHamma(newPointCoordinate.getX(), newPointCoordinate.getY(),
                        oldPoint.getHamma()));
            }
        }
    }

    public double getPhi(Point point) {
        double toReturn = point.getX() * cos(alpha) + point.getY() * sin(alpha);

        double sum = 0;
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            sum += shape.getAllPoints().get(i).getHamma() * atan2(
                    point.getY() - shape.getAllPoints().get(i).getY(),
                    point.getX() - shape.getAllPoints().get(i).getX()
            );
        }

        for (List<PointWithHamma> pointWithHammaList : flyingPointsWithHama.values()) {
            for (PointWithHamma pointWithHamma : pointWithHammaList) {
                sum += pointWithHamma.getHamma() * atan2(
                        point.getY() - pointWithHamma.getY(),
                        point.getX() - pointWithHamma.getX()
                );
            }
        }
        toReturn += 1 / (2 * PI) * sum;
        return toReturn;
    }

    public Pair<Double> computeForFlyingPoints(double partialHama, Point point, int firstFlyingPointNumber,
                                               int numberOfLastPoint) {
        Point firstStrangeValue = new Point(0, 0);
        Point secondStrangeValue = new Point(0, 0);
        double toReturn = 0;
        List<PointWithHamma> flyingPoints = flyingPointsWithHama.get(shape.getPointsOfSeparation().get(firstFlyingPointNumber));

        for (int i = firstFlyingPointNumber; i < flyingPoints.size(); i++) {
            PointWithHamma flyingPoint = flyingPoints.get(i);
            partialHama += flyingPoint.getHamma();

            if (i + 1 < flyingPointsWithHama.size()) {
                Point neighborPoint = flyingPoints.get(i + 1);
                firstStrangeValue = neighborPoint.minus(flyingPoint);
            } else {
                Point neighborPoint = shape.getListOfPoints().get(numberOfLastPoint);
                firstStrangeValue = neighborPoint.minus(flyingPoint);
            }

            secondStrangeValue = point.minus(flyingPoint);

            double strangeR = Math.pow(secondStrangeValue.getX(), 2) + Math.pow(secondStrangeValue.getY(), 2);
            strangeR = strangeR < DELTA ? DELTA : strangeR;

            double strangeProd = (firstStrangeValue.getY() * secondStrangeValue.getX() - firstStrangeValue.getX() * secondStrangeValue.getY())
                    / (strangeR * (2 * PI));
            strangeProd *= partialHama;

            toReturn += strangeProd;
        }

        return new Pair<Double>(toReturn, partialHama);
    }

    public Pair<Double> computeForContourPoints(double partialHama, Point point,
                                                int firstPointNumber, int lastPointNumber) {
        Point firstStrangeValue = new Point(0, 0);
        Point secondStrangeValue = new Point(0, 0);
        double toReturn = 0;

        for (int i = firstPointNumber; i < lastPointNumber; i += 1) {
            Point pointOnContour = shape.getAllPoints().get(i);
            partialHama += shape.getAllPoints().get(i).getHamma();

            Point neighborPoint = shape.getAllPoints().get(i + 1);
            firstStrangeValue = neighborPoint.minus(pointOnContour);

            secondStrangeValue = point.minus(pointOnContour);

            double strangeR = Math.pow(secondStrangeValue.getX(), 2) + Math.pow(secondStrangeValue.getY(), 2);
            strangeR = strangeR < DELTA ? DELTA : strangeR;

            double strangeProd = (firstStrangeValue.getY() * secondStrangeValue.getX() - firstStrangeValue.getX() * secondStrangeValue.getY())
                    / (strangeR * (2 * PI));
            strangeProd *= partialHama;

            toReturn += strangeProd;
        }

        return new Pair<Double>(toReturn, partialHama);
    }

    public double getPhiDipol(Point point) {       //todo refactor to use two methods : compute for
        double toReturn = point.getX() * cos(alpha) + point.getY() * sin(alpha);

        Point firstStrangeValue = null;
        Point secondStrangeValue = null;

        double partialHama = 0.0;
        int amountOfPointsOfSeparation = shape.getPointsOfSeparation().size();

        Pair<Double> result = computeForFlyingPoints(0, point, 0, 0);
        toReturn += result.getFirstValue();
        partialHama += result.getSecondValue();

        result = computeForContourPoints(partialHama, point, 0, shape.getAllPoints().size() / 4 - 1);
        toReturn += result.getFirstValue();
        partialHama += result.getSecondValue();

        result = computeForFlyingPoints(0, point, 1, 1);
        toReturn += result.getFirstValue();
        partialHama += result.getSecondValue();

        int middleIndex = shape.getAllPoints().size() / 4 + (shape.getAllPoints().size() / 4 / 2);

        result = computeForContourPoints(partialHama, point, shape.getAllPoints().size() / 4, middleIndex);
        toReturn += result.getFirstValue();
        partialHama += result.getSecondValue();

        double newHama = 0;
        result = computeForFlyingPoints(0, point, 4, 5);
        toReturn += result.getFirstValue();
        newHama += result.getSecondValue();

        result = computeForContourPoints(newHama, point, shape.getAllPoints().size() / 4 * 3, shape.getAllPoints().size() - 1);
        toReturn += result.getFirstValue();
        newHama += result.getSecondValue();

        partialHama += newHama;

        result = computeForContourPoints(partialHama, point, middleIndex, shape.getAllPoints().size() / 4 * 2);
        toReturn += result.getFirstValue();
        partialHama += result.getSecondValue();
        partialHama = 0;

        result = computeForFlyingPoints(0, point, 2, 2);
        toReturn += result.getFirstValue();
        partialHama += result.getSecondValue();

        result = computeForContourPoints(partialHama, point, shape.getAllPoints().size() / 4 * 2, shape.getAllPoints().size() / 4 * 3);
        toReturn += result.getFirstValue();
        partialHama += result.getSecondValue();

       /* for (int i = flyingPointsWithHama.size() - 2; i > 3; i -= amountOfPointsOfSeparation) {
            PointWithHamma flyingPoint = flyingPointsWithHama.get(i);
            partialHama += flyingPoint.getHamma();

            if (i != flyingPointsWithHama.size() - 2) {
                Point neighborPoint = flyingPointsWithHama.get(i - amountOfPointsOfSeparation);
                firstStrangeValue = flyingPoint.minus(neighborPoint);
            } else {
                Point neighborPoint = shape.getListOfPoints().get(3);
                firstStrangeValue = flyingPoint.minus(neighborPoint);
            }

            secondStrangeValue = point.minus(flyingPoint);

            double strangeR = Math.pow(secondStrangeValue.getX(), 2) + Math.pow(secondStrangeValue.getY(), 2);
            strangeR = strangeR < DELTA ? DELTA : strangeR;

            double strangeProd = (firstStrangeValue.getY() * secondStrangeValue.getX() - firstStrangeValue.getX() * secondStrangeValue.getY())
                    / (strangeR * (2 * PI));
            strangeProd *= partialHama;

            toReturn += strangeProd;
        }           */

        return toReturn;
    }

    public double getPsi(Point point) {
        double toReturn = cos(alpha) * point.getY() - sin(alpha) * point.getX();

        double sum = 0;
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            sum += shape.getAllPoints().get(i).getHamma() * log(computeRj(point, shape.getAllPoints().get(i)));
        }
        for (List<PointWithHamma> flyingPoints : flyingPointsWithHama.values()) {
            for (PointWithHamma flyingPoint : flyingPoints) {
                sum += flyingPoint.getHamma() * log(computeRj(point, flyingPoint));
            }
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
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            sum = sum.add(computeVj(point, shape.getAllPoints().get(i)).multiply(shape.getAllPoints().get(i).getHamma()));
        }
        for (List<PointWithHamma> flyingPoints : flyingPointsWithHama.values()) {
            for (PointWithHamma flyingPoint : flyingPoints) {
                sum = sum.add(computeVj(point, flyingPoint).multiply(flyingPoint.getHamma()));
            }
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

    /**
     * Reflects new point if it crosses the carcases
     *
     * @param oldPoint
     * @param newPoint
     * @return
     */
    private Point reflectIsCrossCarcases(Point oldPoint, Point newPoint) {
        boolean lResult = true;
        int j = 0;
        while (j < shape.getAllPoints().size() - 1) {
            if (Point.intersects(newPoint, oldPoint, shape.getAllPoints().get(j), shape.getAllPoints().get(j + 1))) {
                newPoint = Point.reflect(shape.getAllPoints().get(j), shape.getAllPoints().get(j + 1), newPoint);
                j = 0;
                lResult = false;
            } else j++;
        }
        return newPoint;
    }

    public Point awayFromSegmentationPoint(Point point) {
        boolean lResult = false;
        double lAbsoluteDistance = shape.getDelta();
        double lDistance = 0;
        double lDistance2 = lAbsoluteDistance + shape.getDelta();
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            lDistance = 0;
            lDistance2 = lAbsoluteDistance + shape.getDelta();

            lDistance = Point.distance(shape.getAllPoints().get(i), point);
            if (i != shape.getAllPoints().size() - 1)
                lDistance2 = Point.distance(shape.getAllPoints().get(i + 1), point);
            if (lDistance < lAbsoluteDistance) {
                lResult = true;

                if (lDistance2 < lAbsoluteDistance) {
                    Point.StraightLine lLine =
                            new Point.StraightLine(shape.getAllPoints().get(i), shape.getAllPoints().get(i + 1));
                    point = lLine.awayFromLine(point, lAbsoluteDistance);
                } else
                    point = (point.minus(shape.getAllPoints().get(i))).multiply(lAbsoluteDistance / lDistance).add(shape.getAllPoints().get(i));
                //i = 0;
                //break;
            }

        }
        return point;
    }

    private static class Pair<T> {
        T firstValue;
        T secondValue;

        private Pair(T firstValue, T secondValue) {
            this.firstValue = firstValue;
            this.secondValue = secondValue;
        }

        private T getFirstValue() {
            return firstValue;
        }

        private void setFirstValue(T firstValue) {
            this.firstValue = firstValue;
        }

        private T getSecondValue() {
            return secondValue;
        }

        private void setSecondValue(T secondValue) {
            this.secondValue = secondValue;
        }
    }
}
