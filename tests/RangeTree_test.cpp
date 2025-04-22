/***
 * MIT License
 *
 * Copyright (c) 2016 Luca Weihs
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <vector>
#include "RangeTree.h"
#include <random>
#include <cmath>
#include <doctest.h>

namespace RT = Rangetree;

template <typename Scalar>
std::vector<Scalar> abs(std::vector<Scalar> vec) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i] = std::abs(vec[i]);
    }
    return vec;
}

template <typename Scalar>
std::vector<Scalar> add(const std::vector<Scalar>& v0, const std::vector<Scalar>& v1) {
    std::vector<Scalar> v;
    for (int i = 0; i < v0.size(); i++) {
        v.push_back(v0[i] + v1[i]);
    }
    return v;
}

template <typename Scalar, class Value>
std::vector<RT::Point<Scalar,Value>> slice(std::vector<RT::Point<Scalar,Value>> points, int start, int end) {
    std::vector<RT::Point<Scalar,Value>> newVec(points.begin() + start, points.begin() + end + 1);
    return newVec;
}

template <typename Scalar, class Value>
std::vector<RT::Point<Scalar,Value>> sortPoints(std::vector<RT::Point<Scalar,Value>> points) {
    RT::PointOrdering<Scalar,Value> pointOrdering(0);
    std::sort(points.begin(), points.end(), pointOrdering);
    return points;
}

template <typename Scalar, class Value>
std::vector<RT::Point<Scalar,Value>> sortAndMerge(std::vector<RT::Point<Scalar,Value>> points) {
    if (points.size() == 0) {
        return points;
    }
    RT::PointOrdering<Scalar,Value> pointOrdering(0);
    points = sortPoints(points);
    std::vector<RT::Point<Scalar,Value>> sortedPoints = {points[0]};
    int k = 0;
    for (int i = 1; i < points.size(); i++) {
        if (pointOrdering.equals(sortedPoints[k], points[i])) {
            sortedPoints[k].increaseCountBy(points[i].count());
        } else {
            sortedPoints.push_back(points[i]);
            k++;
        }
    }
    return sortedPoints;
}

template <typename Scalar, class Value>
void printPoints(std::vector<RT::Point<Scalar,Value>> points) {
    for (int i = 0; i < points.size(); i++) {
        points[i].print();
    }
}

TEST_CASE("range_tree_count_test")
{
	SUBCASE("returns_correct_for_one_dim")
	{
		std::vector<double> values = { 1.0,1.0,2.0,2.0,3.0,3.0,3.0,4.0,4.0,5.0 };
		std::vector<RT::Point<double, int>> points = {};

		auto f = [](double a) { std::vector<double> b = { a }; return b; };
		for (int i = 0; i < values.size(); i++) {
			RT::Point<double, int> a(f(values[i]), 0);
			points.push_back(a);
		}

		RT::RangeTree<double, int> rtree(points);

		auto g = [](bool a) { std::vector<bool> b = { a }; return b; };
		CHECK_EQ(rtree.countInRange(f(-12.0), f(30.0), g(true), g(true)), 10);
		CHECK_EQ(rtree.countInRange(f(2.0), f(4.0), g(true), g(true)), 7);
		CHECK_EQ(rtree.countInRange(f(2.0), f(4.0), g(false), g(true)), 5);
		CHECK_EQ(rtree.countInRange(f(2.0), f(4.0), g(true), g(false)), 5);
		CHECK_EQ(rtree.countInRange(f(2.0), f(4.0), g(false), g(false)), 3);
		CHECK_EQ(rtree.countInRange(f(3.0), f(3.0), g(false), g(false)), 0);
		CHECK_EQ(rtree.countInRange(f(3.0), f(4.0), g(false), g(false)), 0);
	}

	SUBCASE("returns_correct_for_two_dim")
	{
		std::vector<double> v1 = { 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0, 3.1 };
		std::vector<double> v2 = { 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 3.2 };
		std::vector<RT::Point<double, int>> points = {};

		auto f = [](double a, double b) { std::vector<double> c = { a, b }; return c; };
		for (int i = 0; i < v1.size(); i++) {
			RT::Point<double, int> a(f(v1[i], v2[i]), 0);
			points.push_back(a);
		}

		RT::RangeTree<double, int> rtree(points);

		auto g = [](bool a, bool b) { std::vector<bool> c = { a, b }; return c; };
		// Selecting 2-dim regions with boundary
		CHECK_EQ(rtree.countInRange(f(0.0, 0.0), f(4.0, 4.0),
			g(true, true), g(true, true)), 11);
		CHECK_EQ(rtree.countInRange(f(0.0, 0.0), f(1.0, 1.0),
			g(true, true), g(true, true)), 2);
		CHECK_EQ(rtree.countInRange(f(0.0, 0.0), f(1.0, 3.0),
			g(true, true), g(true, true)), 4);
		CHECK_EQ(rtree.countInRange(f(1.9, 1.9), f(2.1, 2.1),
			g(true, true), g(true, true)), 1);
		CHECK_EQ(rtree.countInRange(f(-1000.0, -1000.0), f(1000.0, 0.9999999),
			g(true, true), g(true, true)), 0);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(true, true), g(true, true)), 10);

		// Selecting 2-dim regions without (some) boundaries
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(false, true), g(true, true)), 6);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(true, false), g(true, true)), 6);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(true, true), g(false, true)), 7);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(true, true), g(true, false)), 7);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(false, false), g(true, true)), 4);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(true, true), g(false, false)), 5);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(true, false), g(true, false)), 3);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(false, true), g(false, true)), 3);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
			g(false, false), g(false, false)), 1);

		// Selecting 0/1-dim regions with boundary
		CHECK_EQ(rtree.countInRange(f(1.0, 0.0), f(1.0, 2.0),
			g(true, true), g(true, true)), 3);
		CHECK_EQ(rtree.countInRange(f(0.0, 2.0), f(3.0, 2.0),
			g(true, true), g(true, true)), 3);
		CHECK_EQ(rtree.countInRange(f(0.0, 2.0), f(3.0, 2.0),
			g(true, true), g(true, true)), 3);
		CHECK_EQ(rtree.countInRange(f(1.3, 3.0), f(2.3, 3.0),
			g(true, true), g(true, true)), 1);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
			g(true, true), g(true, true)), 2);
		CHECK_EQ(rtree.countInRange(f(0.0, 0.0), f(0.0, 1.0),
			g(true, true), g(true, true)), 0);

		// Selecting 0/1-dim regions without (some) boundaries
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
			g(false, true), g(true, true)), 0);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
			g(true, false), g(true, true)), 1);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
			g(true, true), g(false, true)), 0);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
			g(true, true), g(true, false)), 2);

		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
			g(false, true), g(true, true)), 0);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
			g(true, false), g(true, true)), 0);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
			g(true, true), g(false, true)), 0);
		CHECK_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
			g(true, true), g(true, false)), 0);
	}

	SUBCASE("returns_correct_for_three_dim")
	{
		auto f = [](double a, double b, double c) { std::vector<double> d = { a, b, c }; return d; };
		std::vector<RT::Point<double, int>> points = {};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					RT::Point<double, int> a(f(i, j, k), 0);
					points.push_back(a);
					if (i == 1 && j == 1 && k == 1) {
						points.push_back(a);
					}
				}
			}
		}
		RT::RangeTree<double, int> rtree(points);
		RT::NaiveRangeCounter<double, int> nrc(points);

		auto g = [](bool a, bool b, bool c) { std::vector<bool> d = { a, b, c }; return d; };
		std::vector<std::vector<bool>> boolVectors = {};
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					std::vector<bool> a = {};
					a.push_back(i == 1);
					a.push_back(j == 1);
					a.push_back(k == 1);
					boolVectors.push_back(a);
				}
			}
		}

		// Selecting 0,1,2,3 dim region with all possible boundaries
		auto lower3d = { 0.0, 0.0, 0.0 };
		auto upper3d = { 2.0, 2.0, 2.0 };

		auto lower2d = { 1.0, 0.0, 0.0 };
		auto upper2d = { 1.0, 2.0, 2.0 };

		auto lower1d = { 1.0, 1.0, 0.0 };
		auto upper1d = { 1.0, 1.0, 2.0 };

		auto lower0d = { 1.0, 1.0, 1.0 };
		auto upper0d = { 1.0, 1.0, 1.0 };

		for (int i = 0; i < boolVectors.size(); i++) {
			for (int j = 0; j < boolVectors.size(); j++) {
				auto withLower = boolVectors[i];
				auto withUpper = boolVectors[j];

				CHECK_EQ(nrc.countInRange(lower3d, upper3d, withLower, withUpper),
					rtree.countInRange(lower3d, upper3d, withLower, withUpper));
				CHECK_EQ(nrc.countInRange(lower2d, upper2d, withLower, withUpper),
					rtree.countInRange(lower2d, upper2d, withLower, withUpper));
				CHECK_EQ(nrc.countInRange(lower1d, upper1d, withLower, withUpper),
					rtree.countInRange(lower1d, upper1d, withLower, withUpper));
				CHECK_EQ(nrc.countInRange(lower0d, upper0d, withLower, withUpper),
					rtree.countInRange(lower0d, upper0d, withLower, withUpper));
			}
		}
	}
}

TEST_CASE("range_tree_return_points_test")
{
	SUBCASE("returns_correct_for_one_dim")
	{
		std::vector<double> values = { 3.0, 1.0, 2.0, 11.0, 5.0, 11.0 };
		std::vector<int> counts = { 1, 3, 4, 1, 2, 1 };
		std::vector<double> sortedValues = { 1.0, 2.0, 3.0, 5.0, 11.0 };
		std::vector<int> sortedCounts = { 3, 4, 1, 2, 2 };
		std::vector<RT::Point<double, int>> points = {};
		std::vector<RT::Point<double, int>> sortedPoints = {};

		auto f = [](double a) { std::vector<double> b = { a }; return b; };
		for (int i = 0; i < values.size(); i++) {
			RT::Point<double, int> a(f(values[i]), 0);
			a.increaseCountBy(counts[i] - 1);
			points.push_back(a);
		}

		std::cout << "Sorted points:" << std::endl;
		auto temp = sortAndMerge(points);
		for (int i = 0; i < sortedValues.size(); i++) {
			RT::Point<double, int> a(f(sortedValues[i]), 0);
			a.increaseCountBy(sortedCounts[i] - 1);
			sortedPoints.push_back(a);
		}
		CHECK_EQ(sortedPoints, sortAndMerge(points));

		RT::RangeTree<double, int> rtree(points);

		auto g = [](bool a) { std::vector<bool> b = { a }; return b; };
		CHECK_EQ(sortPoints(rtree.pointsInRange(f(-12.0), f(30.0), g(true), g(true))),
			sortedPoints);
		CHECK_EQ(sortPoints(rtree.pointsInRange(f(2.0), f(3.0), g(true), g(true))),
			slice(sortedPoints, 1, 2));
		CHECK_EQ(sortPoints(rtree.pointsInRange(f(2.0), f(3.0), g(false), g(true))),
			slice(sortedPoints, 2, 2));
		CHECK_EQ(sortPoints(rtree.pointsInRange(f(2.0), f(3.0), g(true), g(false))),
			slice(sortedPoints, 1, 1));
		CHECK_EQ(sortPoints(rtree.pointsInRange(f(2.0), f(3.0), g(false), g(false))),
			slice(sortedPoints, 1, 0));
	}

	SUBCASE("returns_correct_for_three_dim_cube")
	{
		std::vector<RT::Point<double, int>> points = {};
		auto f = [](double a, double b, double c) { std::vector<double> d = { a, b, c }; return d; };
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					RT::Point<double, int> a(f(i, j, k), 0);
					points.push_back(a);
					points.push_back(a);
				}
			}
		}

		RT::RangeTree<double, int> rtree(points);
		RT::NaiveRangeCounter<double, int> nrc(points);

		auto g = [](bool a, bool b) { std::vector<bool> c = { a, b }; return c; };
		std::vector<std::vector<bool>> boolVectors = {};
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					std::vector<bool> a = {};
					a.push_back(i == 1);
					a.push_back(j == 1);
					a.push_back(k == 1);
					boolVectors.push_back(a);
				}
			}
		}

		// Selecting 0,1,2 dim region with all possible boundaries
		auto lower2d = { 1.0, 0.0, 0.0 };
		auto upper2d = { 1.0, 2.0, 2.0 };

		auto lower1d = { 1.0, 1.0, 0.0 };
		auto upper1d = { 1.0, 1.0, 2.0 };

		auto lower0d = { 1.0, 1.0, 1.0 };
		auto upper0d = { 1.0, 1.0, 1.0 };

		for (int i = 0; i < boolVectors.size(); i++) {
			for (int j = 0; j < boolVectors.size(); j++) {
				auto withLower = boolVectors[i];
				auto withUpper = boolVectors[j];

				CHECK_EQ(sortAndMerge(nrc.pointsInRange(lower2d, upper2d, withLower, withUpper)),
					sortPoints(rtree.pointsInRange(lower2d, upper2d, withLower, withUpper)));
				CHECK_EQ(sortAndMerge(nrc.pointsInRange(lower1d, upper1d, withLower, withUpper)),
					sortPoints(rtree.pointsInRange(lower1d, upper1d, withLower, withUpper)));
				CHECK_EQ(sortAndMerge(nrc.pointsInRange(lower0d, upper0d, withLower, withUpper)),
					sortPoints(rtree.pointsInRange(lower0d, upper0d, withLower, withUpper)));
			}
		}
	}

	SUBCASE("returns_correct_for_three_dim")
	{
		auto f = [](double a, double b, double c) { std::vector<double> d = { a, b, c }; return d; };
		std::vector<RT::Point<double, int>> points = {};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					RT::Point<double, int> a(f(i, j, k), 0);
					points.push_back(a);
					if (i == 1 && j == 1 && k == 1) {
						points.push_back(a);
					}
				}
			}
		}
		RT::RangeTree<double, int> rtree(points);
		RT::NaiveRangeCounter<double, int> nrc(points);

		auto g = [](bool a, bool b, bool c) { std::vector<bool> d = { a, b, c }; return d; };
		std::vector<std::vector<bool>> boolVectors = {};
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					std::vector<bool> a = {};
					a.push_back(i == 1);
					a.push_back(j == 1);
					a.push_back(k == 1);
					boolVectors.push_back(a);
				}
			}
		}

		// Selecting 0,1,2,3 dim region with all possible boundaries
		auto lower3d = { 0.0, 0.0, 0.0 };
		auto upper3d = { 2.0, 2.0, 2.0 };

		auto lower2d = { 1.0, 0.0, 0.0 };
		auto upper2d = { 1.0, 2.0, 2.0 };

		auto lower1d = { 1.0, 1.0, 0.0 };
		auto upper1d = { 1.0, 1.0, 2.0 };

		auto lower0d = { 1.0, 1.0, 1.0 };
		auto upper0d = { 1.0, 1.0, 1.0 };

		for (int i = 0; i < boolVectors.size(); i++) {
			for (int j = 0; j < boolVectors.size(); j++) {
				auto withLower = boolVectors[i];
				auto withUpper = boolVectors[j];

				CHECK_EQ(sortAndMerge(nrc.pointsInRange(lower3d, upper3d, withLower, withUpper)),
					sortPoints(rtree.pointsInRange(lower3d, upper3d, withLower, withUpper)));
				CHECK_EQ(sortAndMerge(nrc.pointsInRange(lower2d, upper2d, withLower, withUpper)),
					sortPoints(rtree.pointsInRange(lower2d, upper2d, withLower, withUpper)));
				CHECK_EQ(sortAndMerge(nrc.pointsInRange(lower1d, upper1d, withLower, withUpper)),
					sortPoints(rtree.pointsInRange(lower1d, upper1d, withLower, withUpper)));
				CHECK_EQ(sortAndMerge(nrc.pointsInRange(lower0d, upper0d, withLower, withUpper)),
					sortPoints(rtree.pointsInRange(lower0d, upper0d, withLower, withUpper)));
			}
		}
	}

	SUBCASE("returns_correct_for_high_dim")
	{
		std::random_device rd;
		std::mt19937 gen(rd());

		std::normal_distribution<> norm(0, 1);
		std::bernoulli_distribution bern(0.5);
		std::binomial_distribution<> binom(10, 0.5);

		int numPoints = 2000;
		int dim = 4;
		std::vector<RangeTree::Point<double, int> > normPoints;
		std::vector<RangeTree::Point<int, int> > binomPoints;
		for (int i = 0; i < numPoints; i++) {
			std::vector<double> dPt;
			std::vector<int> iPt;
			for (int j = 0; j < dim; j++) {
				dPt.push_back(norm(gen));
				iPt.push_back(binom(gen));
			}
			normPoints.push_back(RangeTree::Point<double, int>(dPt, 0));
			binomPoints.push_back(RangeTree::Point<int, int>(iPt, 0));
		}
		RangeTree::RangeTree<double, int> normRt(normPoints);
		RangeTree::RangeTree<int, int> binomRt(binomPoints);

		RangeTree::NaiveRangeCounter<double, int> normNrc(normPoints);
		RangeTree::NaiveRangeCounter<int, int> binomNrc(binomPoints);

		int numTrials = 1000;
		std::uniform_int_distribution<> sampler(0, numPoints - 1);
		for (int i = 0; i < numTrials; i++) {
			std::vector<bool> withLower;
			std::vector<bool> withUpper;

			for (int j = 0; j < dim; j++) {
				withLower.push_back(bern(gen) == 1);
				withUpper.push_back(bern(gen) == 1);
			}

			int ind0 = sampler(gen);
			int ind1 = sampler(gen);

			std::vector<double> normLower = normPoints[ind0].asVector();
			std::vector<double> normUpper = add(normLower, abs(normPoints[ind1].asVector()));

			CHECK_EQ(normRt.countInRange(normLower, normUpper, withLower, withUpper),
				normNrc.countInRange(normLower, normUpper, withLower, withUpper));
			CHECK_EQ(sortPoints(normRt.pointsInRange(normLower, normUpper, withLower, withUpper)),
				sortAndMerge(normNrc.pointsInRange(normLower, normUpper, withLower, withUpper)));

			std::vector<int> binomLower = binomPoints[ind0].asVector();
			std::vector<int> binomUpper = add(binomLower, abs(binomPoints[ind1].asVector()));
			CHECK_EQ(binomRt.countInRange(binomLower, binomUpper, withLower, withUpper),
				binomNrc.countInRange(binomLower, binomUpper, withLower, withUpper));
			CHECK_EQ(sortPoints(binomRt.pointsInRange(binomLower, binomUpper, withLower, withUpper)),
				sortAndMerge(binomNrc.pointsInRange(binomLower, binomUpper, withLower, withUpper)));
		}
	}
}
