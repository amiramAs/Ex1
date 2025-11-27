
/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            ans = new double[lx];

            double denom = (xx[0] - xx[1]) * (xx[0] - xx[2]) * (xx[1] - xx[2]);
            double A     = (xx[2] * (yy[1] - yy[0]) + xx[1] * (yy[0] - yy[2]) + xx[0] * (yy[2] - yy[1])) / denom;
            double B     = (xx[2]*xx[2] * (yy[0] - yy[1]) + xx[1]*xx[1] * (yy[2] - yy[0]) + xx[0]*xx[0] * (yy[1] - yy[2])) / denom;
            double C     = (xx[1] * xx[2] * (xx[1] - xx[2]) * yy[0] + xx[2] * xx[0] * (xx[2] - xx[0]) * yy[1] + xx[0] * xx[1] * (xx[0] - xx[1]) * yy[2]) / denom;

            for(int i=1;i<lx;i++) {
                if (i==0){ans[0]=C;}
                if(i==1){ans[1]=B;}
                if(i==2){ans[2]=A;}
            }
        }
		return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
            for(int i=0;i<=Math.max(p1.length,p2.length);i++) {
                double x1=f(p1,i);
                double x2=f(p2,i);
                if(Math.abs(x1-x2)>EPS) {ans=false;}
            }

		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            for(int i=poly.length-1; i>=0; i--) {
                if(poly[i]!=0) {
                    if(i>1) {
                        if (poly[i] > 0) {
                            ans = ans+" +" + poly[i] + "x^" + i;
                        }
                        if (poly[i] <0){
                            ans = ans +" "+ poly[i] + "x^" + i;
                        }
                    }
                    if(i==1){
                        if (poly[i] > 0) {
                            ans = ans+" +" + poly[i] + "x ";
                        }
                        if (poly[i] <0){
                            ans = ans +" "+ poly[i] + "x ";
                        }
                    }
                    if(i==0) {
                        if (poly[i] > 0) {
                            ans = ans+" +" + poly[i] ;
                        }
                        if (poly[i] <0){
                            ans = ans +" "+ poly[i] ;
                        }                    }
                }
            }
		}
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = x1;

        double[] difference = difference(p1, p2);

        ans= root_rec(difference,x1,x2,eps);

		return ans;
	}
	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = 0;

        double step = Math.abs(x1-x2)/numberOfSegments;

        for(double i=x1;i<x2;i=i+step) {
            double j = i+step;
            double y1=f(p,i);
            double y2=f(p,j);

            double lengthX = Math.pow ((i-j),2);
            double lengthy = Math.pow ((y1-y2),2);

            double length=Math.sqrt(lengthX+lengthy);

            ans = ans+length;
        }

		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;

        if(numberOfTrapezoid<2) {
            numberOfTrapezoid=2;
        }

        double step = Math.abs(x1-x2)/numberOfTrapezoid;

        for(double i=x1;i<x2;i=i+step) {

            double j = i+step;

            if(f(p1,i)>f(p2,i) && f(p1,j)<f(p2,j) || f(p1,i)<f(p2,i) && f(p1,j)>f(p2,j)) {
                double same = sameValue(p1,p2,i,j,EPS);

                double height1 = Math.abs(f(p1,i)-f(p2,i));
                double width1 = same-i;

                ans += (height1 * width1) /2;

                double height2 = Math.abs(f(p1,j)-f(p2,j));
                double width2 = j-same;

                ans+=(height2*width2)/2;

            }else {
                double height1 = Math.abs(f(p1, i) - f(p2, i));
                double height2 = Math.abs(f(p1, j) - f(p2, j));

                ans += (height1 + height2) / 2 * step;
            }
        }

        return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0

        String regex = "[ ]";
        String[] poly = p.split(regex);

        int count=0;

        for(int i=0;i<poly.length;i++) {
            if(poly[i]!="") {
                count++;
            }
        }

        ans=new double[count];

        count=0;

        for(int i=0;i<poly.length;i++) {
            if(poly[i]!="") {
                String[] split = poly[i].split("x");
                if (split[0].charAt(0) == '+') {
                    String a = split[0].substring(1);
                    double x = Double.parseDouble(a);
                    ans[ans.length - count-1] = x;
                }

                if (split[0].charAt(0) == '-') {
                    String a = split[0].substring(1);
                    double x = Double.parseDouble(a) * (-1);
                    ans[ans.length - count-1] = x;
                }

                if (split[0].charAt(0) != '-' && split[0].charAt(0) != '+') {
                    String a = split[0];
                    double x = Double.parseDouble(a) ;
                    ans[ans.length - count-1] = x;
                }

                count++;
            }
        }

            return ans;
	}

	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;
        ans=new double [Math.max(p1.length,p2.length)];
        for (int i=0;i<ans.length;i++) {
            if(i<p1.length&&i<p2.length) {
                ans[i] = p1[i] + p2[i];
            }
            if(i>=p1.length) {
                ans[i] = p2[i];
            }
            if(i>=p2.length) {
                ans[i] = p1[i];
            }
        }
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;

        if(p1.length!=0&&p2.length!=0) {
            ans=new double [p1.length+p2.length];

            for (int i=0;i<p1.length;i++) {
                for(int j=0;j<p2.length;j++) {
                    ans[i+j] = ans[i+j] + p1[i]*p2[j];
                }

            }
        }

        return ans;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;
        if(po.length!=0) {ans=new double [po.length-1];}
        for(int i=1;i<po.length;i++) {
            ans[i-1]=po[i]*i;
        }
		return ans;
	}


    public static double[] difference(double[]p1, double[] p2) {
        double[] inverse= new double[p2.length];
        System.arraycopy(p2, 0, inverse, 0, p2.length);

        for(int i=0;i<inverse.length;i++) {
            inverse[i] = (-1) * p2[i];
        }

        double[] ans = add(p1,inverse);
        return ans;
    }
}
