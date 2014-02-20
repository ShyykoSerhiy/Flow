package shape;

import draw.Point;
import draw.PointWithHamma;

import java.util.List;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 14:18
 */
public abstract class Shape {
    /**
     * List of points which describe shape
     */
    protected List<PointWithHamma> listOfPoints;
    /**
     * All points
     */
    protected List<PointWithHamma> allPoints;
    /**
     * Kolokation points
     */
    protected List<Point> kolokPoints;
    /**
     * Normals in kolocation points
     */
    protected List<Point> normal;
    /**
     * Points from which flying points will take off
     */
    protected List<PointWithHamma> pointsOfSeparation;
    protected double delta;

    public List<PointWithHamma> getListOfPoints() {
        return listOfPoints;
    }

    public List<PointWithHamma> getAllPoints() {
        return allPoints;
    }

    public List<Point> getKolokPoints() {
        return kolokPoints;
    }

    public List<Point> getNormal() {
        return normal;
    }

    public List<PointWithHamma> getPointsOfSeparation() {
        return pointsOfSeparation;
    }

    public double getDelta() {
        return delta;
    }
}
