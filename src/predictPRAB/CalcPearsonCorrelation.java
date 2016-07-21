package predictPRAB;

public class CalcPearsonCorrelation {

	/**
	 * calculate pearson correlation
	 * 
	 * @param xs
	 * @param ys
	 * @return
	 */
	public static double pcc(double[] xs, double[] ys) {
		double sx = 0.0;
		double sy = 0.0;
		double sxx = 0.0;
		double syy = 0.0;
		double sxy = 0.0;

		int n = ys.length;

		for (int i = 0; i < n; i++) {
			double x = xs[i];
			double y = ys[i];

			sx += x;
			sy += y;
			sxx += x * x;
			syy += y * y;
			sxy += x * y;
		}

		double cov = sxy / n - sx * sy / n / n;
		double sigmax = Math.sqrt(sxx / n - sx * sx / n / n);
		double sigmay = Math.sqrt(syy / n - sy * sy / n / n);

		return (cov / sigmax / sigmay);
	}
	
}
