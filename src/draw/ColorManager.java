package draw;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 19:54
 */
public class ColorManager {
    private static final int COLOR_STEPS = 10;

    private double minimalValue;
    private double maximumValue;
    private Color minimumColor;
    private Color maximumColor;

    private double valueStep;
    private Color colorStep;

    public ColorManager(double minimalValue, double maximumValue, Color minimumColor, Color maximumColor) {
        this.minimalValue = minimalValue;
        this.maximumValue = maximumValue;
        this.minimumColor = minimumColor;
        this.maximumColor = maximumColor;

        valueStep = (maximumValue - minimalValue) / COLOR_STEPS;
        colorStep = maximumColor.minus(minimumColor).divide(COLOR_STEPS);
    }

    public Color getColor(double value) {
        double normalizedValue = value - minimalValue;
        int stepNumberOfValue = (int) Math.ceil(normalizedValue / valueStep);

        return minimumColor.plus(colorStep.multiply(stepNumberOfValue));
    }
}
