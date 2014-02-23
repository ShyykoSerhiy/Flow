package shape;

import draw.*;
import draw.Color;
import draw.Point;
import processing.core.PApplet;
import solver.Solver;

import java.util.List;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 14:18
 */
public class ShapeDrawer {
    private static final int POINT_SIZE = 5;
    private static final Color RED = new Color(255, 0, 0);
    private static final Color GREEN = new Color(0, 255, 0);
    private static final Color BLUE = new Color(0, 0, 255);

    private Helper helper;
    private Shape shape;
    private PApplet drawer;
    private Solver solver;

    public ShapeDrawer(Helper helper, Shape shape, PApplet drawer, Solver solver) {
        this.helper = helper;
        this.shape = shape;
        this.drawer = drawer;
        this.solver = solver;

        phiValues = new double[0][0];
        vectorValues = new Point[0][0];

        negativePointColorManager = new ColorManager(0, 10, new Color(0, 0, 0), new Color(0, 0, 0));
        positivePointColorManager = new ColorManager(0, 10, new Color(0, 0, 0), new Color(0, 0, 0));
    }

    public void draw() {
        //drawPoints(shape.getAllPoints(), BLUE);
        //drawPoints(shape.getKolokPoints(), RED);
        //drawPoints(shape.getListOfPoints(), GREEN);
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            drawer.strokeWeight(POINT_SIZE);

            Color color = getColorByHamma(shape.getAllPoints().get(i).getHamma());
            drawer.stroke(color.getR(), color.getG(), color.getB());
            Point drawPoint = helper.getDrawPoint(shape.getAllPoints().get(i));
            drawer.point((float) drawPoint.getX(), (float) drawPoint.getY());
        }
    }

    private void drawPoints(List<? extends Point> points, Color color) {
        drawer.strokeWeight(POINT_SIZE);
        drawer.stroke(color.getR(), color.getG(), color.getB());
        for (Point point : points) {
            drawPoint(point);
        }
    }

    private void drawPoint(Point coordinatePoint) {
        Point drawPoint = helper.getDrawPoint(coordinatePoint);
        drawer.point((float) drawPoint.getX(), (float) drawPoint.getY());
    }


    //todo move to another class
    private double[][] phiValues;
    private ColorManager colorManager;
    private static final int POINT_STEP = 5;

    public void computePhi() {
        phiValues = new double[helper.getWidth()][helper.getHeight()];

        double maxPhi = Double.NEGATIVE_INFINITY;
        double minPhi = Double.POSITIVE_INFINITY;
        for (int i = 0; i < phiValues.length; i += POINT_STEP) {
            for (int j = 0; j < phiValues[i].length; j += POINT_STEP) {
                Point point = helper.getCoordinatePoint(new Point(i, j));
                double phi = solver.getPsi(point);
                if (phi < minPhi) {
                    minPhi = phi;
                }
                if (phi > maxPhi) {
                    maxPhi = phi;
                }
                phiValues[i][j] = phi;
            }
        }
        colorManager = new ColorManager(minPhi, maxPhi, new Color(152, 255, 152), new Color(178, 34, 34));
    }

    public void drawPhi() {
        drawer.strokeWeight(POINT_STEP * 2);
        for (int i = 0; i < phiValues.length; i += POINT_STEP) {
            for (int j = 0; j < phiValues[i].length; j += POINT_STEP) {
                Color color = colorManager.getColor(phiValues[i][j]);
                drawer.stroke(color.getR(), color.getG(), color.getB());
                drawer.point((float) i, (float) j);
            }
        }
    }

    private Point[][] vectorValues;

    public void computeVectors() {
        vectorValues = new Point[helper.getWidth()][helper.getHeight()];
        for (int i = 0; i < vectorValues.length; i += 10) {
            for (int j = 0; j < vectorValues[i].length; j += 10) {
                vectorValues[i][j] = solver.computeV(helper.getCoordinatePoint(new Point(i, j))).multiply(5);
            }
        }
    }

    public void drawVectors() {
        drawer.stroke(10, 10, 10);
        drawer.strokeWeight(2);

        for (int i = 0; i < vectorValues.length; i += 10) {
            for (int j = 0; j < vectorValues[i].length; j += 10) {
                Point vector = vectorValues[i][j].add(new Point(i, j));
                drawer.line((float) i, (float) j, (float) vector.getX(), (float) vector.getY());
            }
        }
    }

    private ColorManager negativePointColorManager;
    private ColorManager positivePointColorManager;

    public void computeColorManagers() {
        double maximumValue = 0;
        double minimumValue = 0;

        for (List<PointWithHamma> pointWithHammaList : solver.getFlyingPointsWithHama().values()) {
            for (PointWithHamma pointWithHamma : pointWithHammaList) {
                if (pointWithHamma.getHamma() > maximumValue) {
                    maximumValue = pointWithHamma.getHamma();
                }
                if (pointWithHamma.getHamma() < minimumValue) {
                    minimumValue = pointWithHamma.getHamma();
                }
            }
        }

        for (PointWithHamma pointWithHamma : shape.getAllPoints()) {
            double hamma = pointWithHamma.getHamma();
            if (hamma > maximumValue) {
                maximumValue = hamma;
            }
            if (hamma < minimumValue) {
                minimumValue = hamma;
            }
        }

        negativePointColorManager = new ColorManager(minimumValue, 0, new Color(10, 10, 255), new Color(255, 255, 255));
        positivePointColorManager = new ColorManager(0, maximumValue, new Color(255, 255, 255), new Color(255, 10, 10));
    }

    public void drawFlyingPoints() {
        for (List<PointWithHamma> pointWithHammaList : solver.getFlyingPointsWithHama().values()) {
            for (PointWithHamma pointWithHamma : pointWithHammaList) {
                drawer.strokeWeight(POINT_SIZE);

                Color color = getColorByHamma(pointWithHamma.getHamma());
                drawer.stroke(color.getR(), color.getG(), color.getB());
                Point drawPoint = helper.getDrawPoint(pointWithHamma);
                drawer.point((float) drawPoint.getX(), (float) drawPoint.getY());
            }
        }
    }

    public Color getColorByHamma(double hamma) {
        Color color;
        if (hamma < 0) {
            color = negativePointColorManager.getColor(hamma);
        } else {
            color = positivePointColorManager.getColor(hamma);
        }
        return color;
    }
}
