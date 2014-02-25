package shape;

import draw.Point;
import draw.PointWithHamma;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * User: shyyko
 * Date: 02.02.14
 * Time: 19:18
 */
public class ShShape extends Shape {
    private static final int AMOUNT_OF_SIDES = 4;
    private int amountOfAdditionalPoints;

    public ShShape(int amountOfAdditionalPoints) {
        assert amountOfAdditionalPoints % AMOUNT_OF_SIDES == 0;
        this.amountOfAdditionalPoints = amountOfAdditionalPoints;

        listOfPoints = new ArrayList<PointWithHamma>(6);
        listOfPoints.addAll(Arrays.asList(
                new PointWithHamma(2.0f, 1.0f, 0), new PointWithHamma(2.0f, 2.0f, 0), new PointWithHamma(1.0f, 2.0f, 0), new PointWithHamma(1.0f, 1.0f, 0),
                new PointWithHamma(1.5f, 1.0f, 0), new PointWithHamma(1.5f, 2.0f, 0)
        ));
        pointsOfSeparation = new ArrayList<PointWithHamma>(5);
        pointsOfSeparation.addAll(Arrays.asList(
                listOfPoints.get(0), listOfPoints.get(1), listOfPoints.get(2), listOfPoints.get(3),
                listOfPoints.get(4)
        ));

        kolokPoints = new ArrayList<Point>();

        allPoints = new ArrayList<PointWithHamma>();
        allPoints.add(listOfPoints.get(0));
        allPoints.addAll(createAdditionalPoints(listOfPoints.get(0), listOfPoints.get(1), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        allPoints.add(listOfPoints.get(1));
        allPoints.addAll(createAdditionalPoints(listOfPoints.get(1), listOfPoints.get(2), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        //we need to add one more point (middle point to bind middle contour)
        allPoints.add(listOfPoints.get(2));
        allPoints.addAll(createAdditionalPoints(listOfPoints.get(2), listOfPoints.get(3), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        allPoints.add(listOfPoints.get(3));

        kolokPoints.addAll(createCollocPoints(allPoints));
        //additional side(middle)
        List<PointWithHamma> middlePoints = new ArrayList<PointWithHamma>(amountOfAdditionalPoints / AMOUNT_OF_SIDES + 2);
        middlePoints.add(listOfPoints.get(4));
        middlePoints.addAll(createAdditionalPoints(listOfPoints.get(4), listOfPoints.get(5), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        middlePoints.add(listOfPoints.get(5));
        // we need to replace the nearest point by listOfPoints.get(5), so it will have the same hamma
        /*double minDistance = Integer.MAX_VALUE;
        int minIndex = 0;
        PointWithHamma middlePoint = listOfPoints.get(5);
        for (int i = 0; i < allPoints.size(); i++) {
            double distance = Point.distance(allPoints.get(i), middlePoint);
            if (distance < minDistance){
                minDistance = distance;
                minIndex = i;
            }
        }
        allPoints.set(minIndex, middlePoint);*/
        allPoints.addAll(middlePoints);

        kolokPoints.addAll(createCollocPoints(middlePoints));


        normal = new ArrayList<Point>();
        int kolokPointsPerSide = kolokPoints.size() / AMOUNT_OF_SIDES;
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(1, 0)));
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(0, 1)));
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(-1, 0)));
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(1, 0)));

        assert normal.size() == kolokPoints.size();
        delta = allPoints.get(1).minus(allPoints.get(0)).getY();
    }

    private List<PointWithHamma> createAdditionalPoints(Point firstPoint, Point secondPoint, int amountOfPoints) {
        List<PointWithHamma> additionalPoints = new ArrayList<PointWithHamma>(amountOfPoints);

        Point step = secondPoint.minus(firstPoint).divide(amountOfPoints + 1);
        for (int i = 0; i < amountOfPoints; i++) {
            Point pointToAdd = firstPoint.add(step.multiply(i + 1));
            additionalPoints.add(new PointWithHamma(pointToAdd, 0));
        }
        return additionalPoints;
    }

    private List<Point> createCollocPoints(List<? extends Point> points) {
        List<Point> kolokPoints = new ArrayList<Point>();
        for (int i = 1; i < points.size(); i++) {
            Point step = points.get(i).minus(points.get(i - 1)).divide(2);
            kolokPoints.add(points.get(i - 1).add(step));
        }
        return kolokPoints;
    }
}
