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
#include <doctest.h>

namespace RT = Rangetree;
typedef RT::Point<double, int> Point;
typedef RT::Point<double, int>* PtrToPoint;
typedef RT::SortedPointMatrix<double, int> SortedPointMatrix;

TEST_CASE("sorted_points_matrix_test")
{
    SUBCASE("works_with_1_dim")
    {
        std::vector<double> values = { 3.0, 1.0, 2.0, 11.0, 5.0, 11.0 };
        std::vector<int> counts = { 1, 3, 4, 1, 2, 1 };
        std::vector<double> sortedValues = { 1.0, 2.0, 3.0, 5.0, 11.0 };
        std::vector<int> sortedCounts = { 3, 4, 1, 2, 2 };
        std::vector<PtrToPoint> points = {};

        auto f = [](double a) { std::vector<double> b = { a }; return b; };
        for (int i = 0; i < values.size(); i++) {
            Point* a = new Point(f(values[i]), (int)values[i] + 1);
            a->increaseCountBy(counts[i] - 1);
            points.push_back(PtrToPoint(a));
        }

        SortedPointMatrix spm(points);
        auto sortedPoints = spm.getSortedPointsAtCurrentDim();

        CHECK_EQ(spm.numUniquePoints(), sortedValues.size());
        CHECK_EQ(sortedPoints.size(), sortedValues.size());

        for (int i = 0; i < spm.numUniquePoints(); i++) {
            CHECK_EQ(sortedPoints[i]->asVector(), f(sortedValues[i]));
            CHECK_EQ(sortedPoints[i]->value(), sortedValues[i] + 1);
            CHECK_EQ(sortedPoints[i]->count(), sortedCounts[i]);
        }

        CHECK_EQ(spm.getMidPoint()->asVector(), f(3.0));
        auto pairOfSPMs = spm.splitOnMid();
        auto leftSPM = pairOfSPMs.first;
        auto rightSPM = pairOfSPMs.second;
        CHECK_EQ(leftSPM.numUniquePoints(), 3);
        CHECK_EQ(rightSPM.numUniquePoints(), 2);
        CHECK_EQ(leftSPM.getSortedPointsAtCurrentDim()[0]->asVector()[0], 1.0);
        CHECK_EQ(leftSPM.getSortedPointsAtCurrentDim()[1]->asVector()[0], 2.0);
        CHECK_EQ(leftSPM.getSortedPointsAtCurrentDim()[2]->asVector()[0], 3.0);
        CHECK_EQ(leftSPM.getSortedPointsAtCurrentDim()[1]->count(), 4);
        CHECK_EQ(leftSPM.getSortedPointsAtCurrentDim()[1]->value(), 3.0);
    }

    SUBCASE("works_with_3_dim")
    {
        std::vector<double> pt0 = { 3.0, 4.0, 1.0 };
        std::vector<double> pt1 = { 5.0, 2.0, 8.0 };
        std::vector<double> pt2 = { 2.0, 9.0, 3.0 };
        std::vector<double> pt3 = { 3.0, 4.0, 1.0 };
        std::vector<std::vector<double> > pts = { pt0, pt1, pt2, pt3 };

        std::vector<int> values = { 1,2,3,1 };
        std::vector<int> counts = { 1,3,1,1 };

        auto f = [](double a) { std::vector<double> b = { a }; return b; };
        std::vector<PtrToPoint> points0 = {};
        std::vector<PtrToPoint> points1 = {};
        for (int i = 0; i < values.size(); i++) {
            Point* a0 = new Point(pts[i], values[i]);
            Point* a1 = new Point(pts[i], values[i]);
            a0->increaseCountBy(counts[i] - 1);
            a1->increaseCountBy(counts[i] - 1);
            points0.push_back(PtrToPoint(a0));
            points1.push_back(PtrToPoint(a1));
        }

        SortedPointMatrix spm(points0);
        std::vector<double> sorted0 = { 2.0, 3.0, 5.0 };
        std::vector<double> sorted1 = { 2.0, 4.0, 9.0 };
        std::vector<double> sorted2 = { 1.0, 3.0, 8.0 };

        auto g = [](std::vector<PtrToPoint> pts, int dim) {
            std::vector<double> a;
            for (int i = 0; i < pts.size(); i++) { a.push_back(pts[i]->asVector()[dim]); }
            return a;
            };

        CHECK_EQ(spm.numUniquePoints(), 3);
        CHECK_EQ(spm.getSortedPointsAtCurrentDim()[1]->count(), 2);
        CHECK_EQ(g(spm.getSortedPointsAtCurrentDim(), 0), sorted0);

        spm.moveToNextDimension();
        CHECK_EQ(g(spm.getSortedPointsAtCurrentDim(), 1), sorted1);

        spm.moveToNextDimension();
        CHECK_EQ(g(spm.getSortedPointsAtCurrentDim(), 2), sorted2);

        spm = SortedPointMatrix(points1);

        auto pairOfSPMs = spm.splitOnMid();
        auto leftSPM = pairOfSPMs.first;
        auto rightSPM = pairOfSPMs.second;

        CHECK_EQ(leftSPM.numUniquePoints(), 2);
        CHECK_EQ(rightSPM.numUniquePoints(), 1);

        std::vector<double> leftSorted0 = { 2.0, 3.0 };
        std::vector<double> leftSorted1 = { 4.0, 9.0 };
        std::vector<double> leftSorted2 = { 1.0, 3.0 };

        CHECK_EQ(leftSPM.numUniquePoints(), 2);
        CHECK_EQ(leftSPM.getSortedPointsAtCurrentDim()[1]->count(), 2);
        CHECK_EQ(g(leftSPM.getSortedPointsAtCurrentDim(), 0), leftSorted0);

        leftSPM.moveToNextDimension();
        CHECK_EQ(g(leftSPM.getSortedPointsAtCurrentDim(), 1), leftSorted1);

        auto leftLeftSPM = leftSPM.splitOnMid().first;
        std::vector<double> leftSorted11 = { 4.0 };
        std::vector<double> leftSorted21 = { 1.0 };
        CHECK_EQ(leftLeftSPM.getMidPoint()->asVector()[1], 4);
        CHECK_EQ(g(leftLeftSPM.getSortedPointsAtCurrentDim(), 1), leftSorted11);
        leftLeftSPM.moveToNextDimension();
        CHECK_EQ(g(leftLeftSPM.getSortedPointsAtCurrentDim(), 2), leftSorted21);

        std::vector<double> rightSorted0 = { 5.0 };
        std::vector<double> rightSorted1 = { 2.0 };
        std::vector<double> rightSorted2 = { 8.0 };

        CHECK_EQ(rightSPM.numUniquePoints(), 1);
        CHECK_EQ(rightSPM.getSortedPointsAtCurrentDim()[0]->count(), 3);
        CHECK_EQ(g(rightSPM.getSortedPointsAtCurrentDim(), 0), rightSorted0);
        rightSPM.moveToNextDimension();
        CHECK_EQ(g(rightSPM.getSortedPointsAtCurrentDim(), 1), rightSorted1);
        rightSPM.moveToNextDimension();
        CHECK_EQ(g(rightSPM.getSortedPointsAtCurrentDim(), 2), rightSorted2);


        /*CHECK_EQ(spm.numUniquePoints(), sortedValues.size());
        CHECK_EQ(sortedPoints.size(), sortedValues.size());

        for (int i = 0; i < spm.numUniquePoints(); i++) {
            CHECK_EQ(sortedPoints[i]->asVector(), f(sortedValues[i]));
            CHECK_EQ(sortedPoints[i]->value(), sortedValues[i] + 1);
            CHECK_EQ(sortedPoints[i]->count(), sortedCounts[i]);
        }

        CHECK_EQ(spm.getMidPoint()->asVector(), f(3.0));
        auto pairOfSPMs = spm.splitOnMid();
        auto leftSPM = pairOfSPMs.first;
        auto rightSPM = pairOfSPMs.second;
        CHECK_EQ(leftSPM.numUniquePoints(), 3);
        CHECK_EQ(rightSPM.numUniquePoints(), 2);
        CHECK_EQ(leftSPM.getSortedPointsAtDimension(0)[0]->asVector()[0], 1.0);
        CHECK_EQ(leftSPM.getSortedPointsAtDimension(0)[1]->asVector()[0], 2.0);
        CHECK_EQ(leftSPM.getSortedPointsAtDimension(0)[2]->asVector()[0], 3.0);
        CHECK_EQ(leftSPM.getSortedPointsAtDimension(0)[1]->count(), 4);
        CHECK_EQ(leftSPM.getSortedPointsAtDimension(0)[1]->value(), 3.0);*/
    }
}


