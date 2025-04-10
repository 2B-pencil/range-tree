/*
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

/**
* Implements the RangeTree data structure.
*
* See the documentation of the class RangeTree for more details regarding the purpose
* of RangeTree's.
*
* This file contains two main abstractions meant for use:
* 1. The RangeTree class.
* 2. A Point class which captures the idea of a d-dimensional euclidean point.
*/

#ifndef RANGETREE_H
#define RANGETREE_H

#include <vector>
#include <iostream>
#include <sstream>
#include <numeric>
#include <type_traits>
#include <deque>
#include <cmath>
#include <algorithm>
#include <memory>

namespace RangeTree {

    /**
    * A point in euclidean space.
    *
    * A class that represents a multi-dimensional euclidean point
    * with some associated value. We allow for each point to have an
    * associated value so that some more information can be stored with
    * each point. Points can also have a multiplicity/count, this corresponds
    * to having several duplicates of the same point.
    */
    template<typename Scalar = double, typename Value = size_t>
    class Point {
        static_assert(std::is_arithmetic<Scalar>::value, "Type Scalar must be numeric");
    private:
        std::vector<Scalar> vec;
        Value val;
        int multiplicity;

    public:
        /**
        * Constructs an empty point.
        *
        * Creates a point in 0 dimensional euclidean space. This constructor
        * is provided only to make certain edge cases easier to handle.
        */
        Point() : multiplicity(0) {}

        /**
        * Constructs a point.
        *
        * Creates a point with its position in euclidean space defined by vec,
        * value defined by val, and a multiplicity/count of 1.
        *
        * @param vec the position in euclidean space.
        * @param val the value associated with the point.
        */
        Point(const std::vector<Scalar>& vec, const Value& val): val(val), vec(vec), multiplicity(1) {}

        /**
        * Constructs a point.
        *
        * Copies a point.
        *
        * @param vec the position in euclidean space.
        * @param val the value associated with the point.
        */
        Point(const Point<Scalar,Value>& p): val(p.val), vec(p.vec), multiplicity(p.count()) {}


        /**
        * Euclidean position of the point.
        *
        * @return the euclidean position of the point as a std::vector.
        */
        const std::vector<Scalar>& asVector() const { return vec; }

        /**
        * The point's ambient dimension.
        *
        * @return the dimension of the space in which the point lives. I.e. a point of the
        *         form (1,2,3) lives in dimension 3.
        */
        inline int dim() const { return static_cast<int>(vec.size()); }

        /**
        * The point's count/multiplicity.
        *
        * @return returns the count/multiplicity.
        */
        inline int count() const {  return multiplicity; }

        /**
        * Increase the point's count/multiplicity.
        *
        * @param n amount to increase by.
        */
        inline void increaseCountBy(const int& n) {
            if (n < 0) {
                throw std::logic_error("Can't increase by a negative amount");
            }
            multiplicity += n;
        }

        /**
        * The point's value.
        *
        * @return the value stored in the point.
        */
        inline Value value() const { return val; }

        /**
        * Index a point.
        *
        * Get the ith coordinate value of the point. I.e. if a point is of the form (4,5,6),
        * then its 0th coordinate value is 4 while its 2nd is 6.
        *
        * @param index the coordinate to index.
        * @return the coordinate value.
        */
        Scalar operator[](int index) const {
            if(index < 0 || index >= static_cast<int>(dim())) {
                throw std::out_of_range("[] access index for point is out of range.");
            }
            return vec[index];
        }

        /**
        * Check for equality.
        *
        * Two points are considered equal if they are in the same spot, have the same
        * multiplicity/count, and store the same value.
        *
        * @param p some other point
        * @return true if \p equals the current point, otherwise false.
        */
        bool operator==(const Point<Scalar,Value>& p) const {
            return vec == p.vec && multiplicity == p.multiplicity && val == p.val;
        }

        /**
        * Check for inequality.
        *
        * The opposite of ==.
        *
        * @param p some other point.
        * @return false if \p equals the current point, otherwise true.
        */
        bool operator!=(const Point<Scalar,Value>& p) const {
            return !((*this) == p);
        }

        /**
        * Prints the point to standard out.
        *
        * As an example, a point with euclidean location (3,4,5) and with a
        * multiplicity/count of 4 will be printed as
        *
        * (3, 4, 5) : 4
        *
        * @param withCount whether or not to display the points count/multiplicity.
        */
        void print(bool withCount=true) const {
            std::cout << "(";
            for (int i = 0; i < dim() - 1; i++) {
                std::cout << (*this)[i] << ", ";
            }
            if (withCount) {
                std::cout << (*this)[dim() - 1] << ") : " << count() << std::endl;
            } else {
                std::cout << (*this)[dim() - 1] << ") : " << std::endl;
            }
        }
    };

    /**
    * A class that totally orders Point<Scalar,Value>'s in euclidean space.
    *
    * A total order of Points is required in the RangeTree. This is an implementation
    * detail that can be ignored. Given a start index \compareStartIndex, this class
    * orders points so that, for two points p_1 = (p_{11}, p_{12},...,p_{1n}) and
    * p_2 = (p_{21}, p_{22},...,p_{2n}) we have that p_1 < p_2 if and only if
    *
    * (p_{1\compareStartInd},...,p_{1n}) < (p_{2\compareStartInd},...,p_{2n})
    *
    * using the usual lexicographic order, or
    *
    * (p_{1\compareStartInd},...,p_{1n}) == (p_{2\compareStartInd},...,p_{2n}) and
    * (p_{11},...,p_{1(\compareStartInd-1)}) < (p_{21},...,p_{2(\compareStartInd-1)})
    *
    * again using the usual lexicographic order.
    */
    template<typename Scalar = double, typename Value = size_t>
    class PointOrdering {
        static_assert(std::is_arithmetic<Scalar>::value, "Type Scalar must be numeric");
    private:
        int compareStartIndex;

    public:
        PointOrdering(int compareStartIndex): compareStartIndex(compareStartIndex) {
            if (compareStartIndex < 0) {
                throw new std::logic_error("Cannot have comparison start index <0.");
            }
        }

        static bool equals(const Point<Scalar,Value>& p1, const Point<Scalar,Value>& p2) {
            return p1.asVector() == p2.asVector();
        }

        inline int getCompareStartIndex() const { return compareStartIndex; }

        bool less(const Point<Scalar,Value>& p1, const Point<Scalar,Value>& p2) const {
            if (p1.dim() != p2.dim()) {
                throw std::logic_error("Points are incomparable (differing dims).");
            }
            if (compareStartIndex >= p1.dim()) {
                throw std::logic_error("Cannot compare points, compare start index >= point dimension.");
            }
            for (int i = compareStartIndex; i < p1.dim(); ++i) {
                if (p1[i] < p2[i]) {
                    return true;
                } else if (p1[i] > p2[i]) {
                    return false;
                }
            }
            for (int i = 0; i < compareStartIndex; ++i) {
                if (p1[i] < p2[i]) {
                    return true;
                } else if (p1[i] > p2[i]) {
                    return false;
                }
            }
            return false;
        }

        bool lessOrEq(const Point<Scalar,Value>& p1, const Point<Scalar,Value>& p2) const {
            return less(p1, p2) || equals(p1, p2);
        }

        bool greater(const Point<Scalar,Value>& p1, const Point<Scalar,Value>& p2) const {
            return less(p2, p1);
        }

        bool greaterOrEq(const Point<Scalar,Value>& p1, const Point<Scalar,Value>& p2) const {
            return greater(p1, p2) || equals(p1,p2);
        }

        bool operator()(const Point<Scalar,Value>& p1, const Point<Scalar,Value>& p2) const {
            return this->less(p1, p2);
        }
    };

    /**
    * A matrix that keeps a collection of points sorted on each coordinate
    */
    template<typename Scalar = double, typename Value = size_t>
    class SortedPointMatrix {
        static_assert(std::is_arithmetic<Scalar>::value, "Type Scalar must be numeric");
    private:
        std::vector<Point<Scalar,Value>* > pointsSortedByCurrentDim;
        std::deque<std::vector<int> > redirectionTable;
        int currentDim;
        int dim;
        static const int MAX_POINTS_BEFORE_SWITCH = 1000;

        std::vector<int> sortOrder(const std::vector<Point<Scalar,Value>* >& points, int onDim) {
            std::vector<int> order(points.size());
            std::iota(order.begin(), order.end(), 0);

            PointOrdering<Scalar,Value> pointOrdering(onDim);
            std::sort(order.begin(), order.end(),
                      [pointOrdering, points](const int& i, const int& j) {
                          return pointOrdering.less(*(points[i]), *(points[j]));
                      });
            return order;
        }

        void sort(std::vector<Point<Scalar, Value>* >& points, const int& onDim) {
            PointOrdering<Scalar,Value> pointOrdering(onDim);
            std::sort(points.begin(), points.end(),
                [pointOrdering](Point<Scalar, Value>* pt0, Point<Scalar, Value>* pt1)
            {
                          return pointOrdering.less(*(pt0), *(pt1));
                      });
        }

        void rearrangeGivenOrder(std::vector<Point<Scalar,Value>* >& points,
                                  const std::vector<int>& order) {
            std::vector<Point<Scalar,Value>* > tmp = points;
            for (size_t i = 0, length = points.size(); i < length; ++i) {
                points[i] = tmp[order[i]];
            }
        }

        SortedPointMatrix(const std::vector<Point<Scalar,Value>* >& pointsSortedByCurrentDim,
                          const std::deque<std::vector<int> >& redirectionTable,
                          int currentDim, int dim) 
            : pointsSortedByCurrentDim(pointsSortedByCurrentDim)
            , redirectionTable(redirectionTable)
            , currentDim(currentDim)
            , dim(dim) {}

    public:
        /**
        * Constructs a sorted point matrix
        */
        SortedPointMatrix(std::vector<Point<Scalar,Value>* >& points)
            : currentDim(0) 
        {
            if (points.size() == 0) {
                throw std::range_error("Cannot construct a SortedPointMatrix with 0 points.");
            }
            else {
                dim = points[0]->dim();
                for (size_t i = 1, length = points.size(); i < length; ++i) {
                    if (points[i]->dim() != dim) {
                        throw std::logic_error("Input points to SortedPointMatrix must all"
                                                       " have the same dimension.");
                    }
                }

                int sortDimension = (points.size() > MAX_POINTS_BEFORE_SWITCH) ? dim - 1 : 0;
                PointOrdering<Scalar,Value> pointOrdering(sortDimension);
                std::sort(points.begin(), points.end(),
                          [pointOrdering](Point<Scalar,Value>* p1, Point<Scalar,Value>* p2) {
                              return pointOrdering.less(*p1, *p2);
                          });
                pointsSortedByCurrentDim.reserve(points.size());
                pointsSortedByCurrentDim.emplace_back(points[0]);
                int k = 0;
                for (size_t i = 1, length = points.size(); i < length; ++i) {
                    if (pointOrdering.equals(*(pointsSortedByCurrentDim[k]), *points[i])) {
                        if (pointsSortedByCurrentDim[k]->value() != points[i]->value()) {
                            throw std::logic_error("Input points have same position but different values");
                        }
                        pointsSortedByCurrentDim[k]->increaseCountBy(points[i]->count());
                    }
                    else {
                        pointsSortedByCurrentDim.emplace_back(points[i]);
                        ++k;
                    }
                }

                if (pointsSortedByCurrentDim.size() > MAX_POINTS_BEFORE_SWITCH) {
                    for (int i = dim - 2; i >= currentDim; --i) {
                        std::vector<int> order = sortOrder(pointsSortedByCurrentDim, i);
                        redirectionTable.emplace_front(order);
                        rearrangeGivenOrder(pointsSortedByCurrentDim, order);
                    }
                }
            }
        }

        void moveToNextDimension() {
            if (currentDim == dim - 1) {
                throw std::logic_error("Already at max dimension, cannot move to next.");
            }
            currentDim++;
            if (pointsSortedByCurrentDim.size() > MAX_POINTS_BEFORE_SWITCH) {
                std::vector<Point<Scalar,Value>* > tmp = pointsSortedByCurrentDim;
                for (size_t i = 0; i < pointsSortedByCurrentDim.size(); ++i) {
                    pointsSortedByCurrentDim[redirectionTable[0][i]] = tmp[i];
                }
                redirectionTable.pop_front();
            } else {
                sort(pointsSortedByCurrentDim, currentDim);
            }
        }

        inline Point<Scalar,Value>* getMidPoint() {
            int mid = (numUniquePoints() - 1) / 2;
            return pointsSortedByCurrentDim[mid];
        }

        inline int numUniquePoints() {
            return  static_cast<int>(pointsSortedByCurrentDim.size());
        }

        inline int getCurrentDim() {
            return currentDim;
        }

        inline std::vector<Point<Scalar,Value>* > getSortedPointsAtCurrentDim() {
            return pointsSortedByCurrentDim;
        }

        /**
        * Constructs two sorted point matrices after splitting on the current midpoint
        */
        std::pair<SortedPointMatrix, SortedPointMatrix> splitOnMid() {
            int n = numUniquePoints();
            if (n == 1) {
                throw std::logic_error("Cannot split on mid when there is only one point.");
            }

            int mid = (n - 1) / 2;
            std::vector<Point<Scalar, Value> *> sortedPointsLeft(mid + 1), sortedPointsRight(n - mid - 1);
            for (int i = 0; i < mid + 1; ++i) {
                sortedPointsLeft[i] = pointsSortedByCurrentDim[i];
            }
            for (int i = mid + 1; i < n; ++i) {
                sortedPointsRight[i - mid - 1] = pointsSortedByCurrentDim[i];
            }

            if (n <= MAX_POINTS_BEFORE_SWITCH) {
                std::deque<std::vector<int> > redirectionTableLeft, redirectionTableRight;
                return std::make_pair(
                        SortedPointMatrix(sortedPointsLeft, redirectionTableLeft, currentDim, dim),
                        SortedPointMatrix(sortedPointsRight, redirectionTableRight, currentDim, dim));
            } else {
                std::vector<bool> onLeft(n);
                for (int i = 0; i < n; i++) {
                    onLeft[i] = i <= mid;
                }

                std::deque<std::vector<int> > redirectionTableLeft(redirectionTable.size(),
                                                                   std::vector<int>(mid + 1));
                std::deque<std::vector<int> > redirectionTableRight(redirectionTable.size(),
                                                                    std::vector<int>(n - mid - 1));
                for (int i = 0, iLength = (int)redirectionTable.size(); i < iLength; ++i) {
                    std::vector<bool> lastOnLeft = onLeft;

                    for (int j = 0, length = numUniquePoints(); j < length; ++j) {
                        onLeft[redirectionTable[i][j]] = lastOnLeft[j];
                    }

                    std::vector<int> newRedirect(numUniquePoints());
                    int kLeft = 0, kRight = 0;
                    for (int j = 0, length = numUniquePoints(); j < length; ++j)
                    {
                        if (onLeft[j])
                        {
                            newRedirect[j] = kLeft;
                            ++kLeft;
                        }
                        else
                        {
                            newRedirect[j] = kRight;
                            ++kRight;
                        }
                    }

                    kLeft = 0, kRight = 0;
                    for (int j = 0, length = numUniquePoints(); j < length; ++j)
                    {
                        if (lastOnLeft[j])
                        {
                            redirectionTableLeft[i][kLeft] = newRedirect[redirectionTable[i][j]];
                            ++kLeft;
                        }
                        else
                        {
                            redirectionTableRight[i][kRight] = newRedirect[redirectionTable[i][j]];
                            ++kRight;
                        }
                    }
                }
                return std::make_pair(
                        SortedPointMatrix(sortedPointsLeft, redirectionTableLeft, currentDim, dim),
                        SortedPointMatrix(sortedPointsRight, redirectionTableRight, currentDim, dim));
            }
        }
    };

    /**
    * A class representing a single node in a RangeTree. These should not be
    * constructed directly, instead use the RangeTree class.
    */
    template<typename Scalar = double, typename Value = size_t>
    class RangeTreeNode {
        static_assert(std::is_arithmetic<Scalar>::value, "Type Scalar must be numeric");
    private:
        std::shared_ptr<RangeTreeNode<Scalar,Value> > left; /**< Contains points <= the comparison point **/
        std::shared_ptr<RangeTreeNode<Scalar,Value> > right; /**< Contains points > the comparison point **/
        std::shared_ptr<RangeTreeNode<Scalar,Value> > treeOnNextDim; /**< Tree on the next dimension **/
        Point<Scalar,Value>* point; /**< The comparison point **/
        bool isLeaf; /**< Whether or not the point is a leaf **/
        int pointCountSum; /**< Total number of points, counting multiplicities, at leaves of the tree **/
        PointOrdering<Scalar,Value> pointOrdering; /**< Helper to totally order input points **/

        // For fractional cascading
        std::vector<Scalar> pointsLastDimSorted;
        std::vector<Point<Scalar,Value>* > allPointsSorted;
        std::vector<int> pointerToGeqLeft;
        std::vector<int> pointerToLeqLeft;
        std::vector<int> pointerToGeqRight;
        std::vector<int> pointerToLeqRight;
        std::vector<int> cumuCountPoints;

		int binarySearchFirstGeq(const Scalar& needle) const {
			int left = 0, right = static_cast<int>(pointsLastDimSorted.size());
			while (left < right) {
				int mid = left + (right - left) / 2;
				if (pointsLastDimSorted[mid] >= needle) {
					right = mid;
				}
				else {
					left = mid + 1;
				}
			}
			return left;
		}

		int binarySearchFirstLeq(const Scalar& needle) const {
			int left = -1, right = static_cast<int>(pointsLastDimSorted.size()) - 1;
			while (left < right) {
				int mid = left + (right - left + 1) / 2;
				if (pointsLastDimSorted[mid] <= needle) {
					left = mid;
				}
				else {
					right = mid - 1;
				}
			}
			return left;
		}

    public:
        /**
        * Construct a range tree structure from points.
        *
        * Creates a range tree structure on the input collection \allPoints using the lexicographic order
        * starting at \compareStartInd.
        *
        * @param uniquePoints a collection of points.
        * @return a range tree structure
        */
        RangeTreeNode(
            SortedPointMatrix<Scalar, Value>& spm,
                      bool onLeftEdge = true,
            bool onRightEdge = true)
            : pointOrdering(spm.getCurrentDim())
        {
            point = spm.getMidPoint();

            if (spm.numUniquePoints() == 1) {
                isLeaf = true;
                pointCountSum = point->count();
                pointsLastDimSorted.emplace_back((*point)[point->dim() - 1]);
                if (spm.getCurrentDim() == point->dim() - 2) {
                    spm.moveToNextDimension();
                }
            }
            else {
                auto spmPair = spm.splitOnMid();
                left = std::shared_ptr<RangeTreeNode<Scalar,Value> >(
                        new RangeTreeNode<Scalar,Value>(spmPair.first, onLeftEdge, false));
                right = std::shared_ptr<RangeTreeNode<Scalar,Value> >(
                        new RangeTreeNode<Scalar,Value>(spmPair.second, false, onRightEdge));
                pointCountSum = left->totalPoints() + right->totalPoints();

                int dim = point->dim();
                if (spm.getCurrentDim() + 2 == dim) {
                    spm.moveToNextDimension();

                    allPointsSorted = spm.getSortedPointsAtCurrentDim();
                    cumuCountPoints.emplace_back(0);
                    for (size_t i = 0, length = allPointsSorted.size(); i < length; ++i) {
                        pointsLastDimSorted.emplace_back((*allPointsSorted[i])[dim - 1]);
                        cumuCountPoints.emplace_back(cumuCountPoints.back() + allPointsSorted[i]->count());
                    }
                    const auto& leftSorted = left->pointsLastDimSorted;
                    const auto& rightSorted = right->pointsLastDimSorted;

                    pointerToGeqLeft = createGeqPointers(pointsLastDimSorted, leftSorted);
                    pointerToGeqRight = createGeqPointers(pointsLastDimSorted, rightSorted);
                    pointerToLeqLeft = createLeqPointers(pointsLastDimSorted, leftSorted);
                    pointerToLeqRight = createLeqPointers(pointsLastDimSorted, rightSorted);
                }
                else if (!onLeftEdge && !onRightEdge && spm.getCurrentDim() + 1 != point->dim()) {
                        spm.moveToNextDimension();
                        treeOnNextDim = std::shared_ptr<RangeTreeNode>(new RangeTreeNode(spm));
                }
                isLeaf = false;
            }
        }

        std::vector<int> createGeqPointers(const std::vector<Scalar>& vec,
                                           const std::vector<Scalar>& subVec) {
            std::vector<int> grePointers(vec.size());
            for (size_t i = 0, length = vec.size(); i < length; ++i) {
                grePointers[i] = static_cast<int>(std::lower_bound(subVec.begin(), subVec.end(), vec[i]) - subVec.begin());
            }
            return grePointers;
        }

        std::vector<int> createLeqPointers(const std::vector<Scalar>& vec,
                                           const std::vector<Scalar>& subVec) {
            std::vector<int> leqPointers(vec.size());
            for (int i = (int)vec.size() - 1; i >= 0; --i) {
				auto it = std::upper_bound(subVec.begin(), subVec.end(), vec[i]);
				leqPointers[i] = static_cast<int>(std::distance(subVec.begin(), it)) - 1;
            }
            return leqPointers;
        }

        /**
        * Construct a RangeTreeNode representing a leaf.
        *
        * @param pointAtLeaf the point to use at the leaf.
        * @param compareStartInd the index defining the lexicographic ordering.
        * @return
        */
        RangeTreeNode(Point<Scalar, Value>* pointAtLeaf, const int& compareStartInd) :
                point(pointAtLeaf), isLeaf(true), pointCountSum(pointAtLeaf->count()), pointOrdering(compareStartInd) {}

        /**
        * Total count of points at the leaves of the range tree rooted at this node.
        *
        * The total count returned INCLUDES the multiplicities/count of the points at the leaves.
        *
        * @return the total count.
        */
        inline int totalPoints() const { return pointCountSum; }

        /**
        * Return all points at the leaves of the range tree rooted at this node.
        * @return all the points.
        */
        std::vector<Point<Scalar,Value> > getAllPoints() const {
            if (isLeaf) {
                std::vector<Point<Scalar,Value> > vec;
                vec.emplace_back(*point);
                return vec;
            }
            auto allPointsLeft = left->getAllPoints();
            auto allPointsRight = right->getAllPoints();

            if(allPointsLeft.size() > allPointsRight.size())
            {
                allPointsLeft.reserve(allPointsLeft.size() + allPointsRight.size());
                allPointsLeft.insert(allPointsLeft.end(), allPointsRight.begin(), allPointsRight.end());
                return allPointsLeft;
            }
            else
            {
                allPointsRight.reserve(allPointsLeft.size() + allPointsRight.size());
                allPointsRight.insert(allPointsRight.end(), allPointsLeft.begin(), allPointsLeft.end());
                return allPointsRight;
            }
        }

        /**
        * Check if point is in a euclidean box.
        *
        * Determines whether or not a point is in a euclidean box. That is, if p_1 = (p_{11},...,p_{1n}) is an
        * n-dimensional point. Then this function returns true if, for all 1 <= i <= n we have
        *
        * lower[i] <= p_{1i} <= upper[i] if withLower[i] == true and withUpper[i] == true, or
        * lower[i] < p_{1i} <= upper[i] if withLower[i] == false and withUpper[i] == true, or
        * lower[i] <= p_{1i} < upper[i] if withLower[i] == true and withUpper[i] == false, or
        * lower[i] < p_{1i} < upper[i] if withLower[i] == false and withUpper[i] == false.
        *
        * @param point the point to check.
        * @param lower the lower points of the rectangle.
        * @param upper the upper bounds of the rectangle.
        * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
        *                  for the lower bounds.
        * @param withUpper as for \withLower but for the upper bounds.
        * @return true if the point is in the rectangle, false otherwise.
        */
        bool pointInRange(const Point<Scalar,Value>& point,
                          const std::vector<Scalar>& lower,
                          const std::vector<Scalar>& upper) const {
            for (int i = 0; i < point.dim(); ++i) {
                if (point[i] < lower[i]) {
                    return false;
                }
                if (point[i] > upper[i]) {
                    return false;
                }
            }
            return true;
        }

        /**
        * Count the number of points at leaves of tree rooted at the current node that are within the given bounds.
        *
        * @param lower see the pointInRange(...) function.
        * @param upper
        * @return the count.
        */
        int countInRange(const std::vector<Scalar>& lower,
                         const std::vector<Scalar>& upper) const {
            if (isLeaf) {
                if (pointInRange(*point, lower, upper)) {
                    return totalPoints();
                }
                else {
                    return 0;
                }
            }
            int compareInd = pointOrdering.getCompareStartIndex();

            if ((*point)[compareInd] > upper[compareInd]) {
                return left->countInRange(lower, upper);
            }
            if ((*point)[compareInd] < lower[compareInd]) {
                return right->countInRange(lower, upper);
            }

            int dim = point->dim();
            if (compareInd + 2 == dim) {
                int n = (int)pointsLastDimSorted.size();
                int geqInd = binarySearchFirstGeq(lower.back());
                int leqInd = binarySearchFirstLeq(upper.back());

                if (geqInd > leqInd) {
                    return 0;
                }
                std::vector<RangeTreeNode<Scalar, Value>* > nodes;
                std::vector<std::pair<int,int> > inds;
                left->leftFractionalCascade(lower,
                                            pointerToGeqLeft[geqInd],
                                            pointerToLeqLeft[leqInd],
                                            nodes,
                                            inds);
                right->rightFractionalCascade(upper,
                                              pointerToGeqRight[geqInd],
                                              pointerToLeqRight[leqInd],
                                              nodes,
                                              inds);
                int sum = 0;
                for (size_t i = 0, length = nodes.size(); i < length; ++i) {
                    if (nodes[i]->isLeaf) {
                        sum += nodes[i]->totalPoints();
                    }
                    else {
                        sum += nodes[i]->cumuCountPoints[inds[i].second + 1] -
                                nodes[i]->cumuCountPoints[inds[i].first];
                    }
                }
                return sum;
            }
            else {
                std::vector<std::shared_ptr<RangeTreeNode<Scalar, Value> > > canonicalNodes;

                if (left->isLeaf) {
                    canonicalNodes.emplace_back(left);
                }
                else {
                    left->leftCanonicalNodes(lower, canonicalNodes);
                }

                if (right->isLeaf) {
                    canonicalNodes.emplace_back(right);
                }
                else {
                    right->rightCanonicalNodes(upper, canonicalNodes);
                }

                int numPointsInRange = 0;
                for (size_t i = 0, length = canonicalNodes.size(); i < length; ++i) {
                    std::shared_ptr<RangeTreeNode<Scalar, Value> > node = canonicalNodes[i];
                    if (node->isLeaf) {
                        if (pointInRange(*(node->point), lower, upper)) {
                            numPointsInRange += node->totalPoints();
                        }
                    }
                    else if (compareInd + 1 == point->dim()) {
                        numPointsInRange += node->totalPoints();
                    }
                    else {
                        numPointsInRange += node->treeOnNextDim->countInRange(lower, upper);
                    }
                }
                return numPointsInRange;
            }
        }

        /**
        * Return the points at leaves of tree rooted at the current node that are within the given bounds.
        *
        * @param lower see the pointInRange(...) function.
        * @param upper
        * @return a std::vector of the Points.
        */
        std::vector<Point<Scalar,Value> > pointsInRange(const std::vector<Scalar>& lower,
                                               const std::vector<Scalar>& upper) const {
            std::vector<Point<Scalar,Value> > pointsToReturn;

            if (isLeaf) {
                if (pointInRange(*point, lower, upper)) {
                    pointsToReturn.emplace_back(*point);
                }
                return pointsToReturn;
            }
            int compareInd = pointOrdering.getCompareStartIndex();

            if ((*point)[compareInd] > upper[compareInd]) {
                return left->pointsInRange(lower, upper);
            }
            if ((*point)[compareInd] < lower[compareInd]) {
                return right->pointsInRange(lower, upper);
            }

            int dim = point->dim();
            if (compareInd + 2 == dim) {
                int n = (int)pointsLastDimSorted.size();
                int geqInd = binarySearchFirstGeq(lower.back());
                int leqInd = binarySearchFirstLeq(upper.back());

                if (geqInd > leqInd) {
                    return pointsToReturn;
                }
                std::vector<RangeTreeNode<Scalar, Value>* > nodes;
                std::vector<std::pair<int,int> > inds;
                left->leftFractionalCascade(lower,
                                            pointerToGeqLeft[geqInd],
                                            pointerToLeqLeft[leqInd],
                                            nodes,
                                            inds);
                right->rightFractionalCascade(upper,
                                              pointerToGeqRight[geqInd],
                                              pointerToLeqRight[leqInd],
                                              nodes,
                                              inds);

                for (size_t i = 0, length = nodes.size(); i < length; ++i)
                {
                    if (nodes[i]->isLeaf)
                    {
                        pointsToReturn.emplace_back(*(nodes[i]->point));
                    }
                    else
                    {
                        for (int j = inds[i].first; j <= inds[i].second; j++)
                        {
                            pointsToReturn.emplace_back(*(nodes[i]->allPointsSorted[j]));
                        }
                    }
                }
                return pointsToReturn;
            }
            else
            {
                std::vector<std::shared_ptr<RangeTreeNode<Scalar, Value> > > canonicalNodes;

                if (left->isLeaf)
                {
                    canonicalNodes.emplace_back(left);
                }
                else
                {
                    left->leftCanonicalNodes(lower, canonicalNodes);
                }

                if (right->isLeaf)
                {
                    canonicalNodes.emplace_back(right);
                }
                else {
                    right->rightCanonicalNodes(upper, canonicalNodes);
                }

                for (int i = 0, length = (int)canonicalNodes.size(); i < length; ++i)
                {
                    std::shared_ptr<RangeTreeNode<Scalar, Value> > node = canonicalNodes[i];
                    if (node->isLeaf)
                    {
                        if (pointInRange(*(node->point), lower, upper))
                        {
                            pointsToReturn.emplace_back(*(node->point));
                        }
                    }
                    else if (compareInd + 1 == point->dim())
                    {
                        auto allPointsAtNode = node->getAllPoints();
                        pointsToReturn.insert(pointsToReturn.end(), allPointsAtNode.begin(), allPointsAtNode.end());
                    }
                    else
                    {
                        auto allPointsAtNode = node->treeOnNextDim->pointsInRange(lower, upper);
                        pointsToReturn.insert(pointsToReturn.end(), allPointsAtNode.begin(), allPointsAtNode.end());
                    }
                }
                return pointsToReturn;
            }
        }

        void leftFractionalCascade(const std::vector<Scalar>& lower,
                                  int geqInd,
                                  int leqInd,
                                  std::vector<RangeTreeNode<Scalar,Value>* >& nodes,
                                  std::vector<std::pair<int,int> >& inds) {
            if (leqInd < geqInd) {
                return;
            }

            int compareInd = point->dim() - 2;

            if (lower[compareInd] <= (*point)[compareInd]) {
                if (isLeaf) {
                    nodes.emplace_back(this);
                    inds.emplace_back(std::pair<int,int>(0,0));
                    return;
                }

                int geqIndRight = pointerToGeqRight[geqInd];
                int leqIndRight = pointerToLeqRight[leqInd];
                if (leqIndRight >= geqIndRight) {
                    nodes.emplace_back(right.get());
                    if (right->isLeaf) {
                        inds.emplace_back(std::pair<int,int>(0,0));
                    }
                    else {
                        inds.emplace_back(std::pair<int,int>(geqIndRight,leqIndRight));
                    }
                }

                left->leftFractionalCascade(lower,
                                            pointerToGeqLeft[geqInd],
                                            pointerToLeqLeft[leqInd],
                                            nodes,
                                            inds);
            }
            else {
                if (isLeaf) {
                    return;
                }
                right->leftFractionalCascade(lower,
                                             pointerToGeqRight[geqInd],
                                             pointerToLeqRight[leqInd],
                                             nodes,
                                             inds);
            }
        }

        void rightFractionalCascade(const std::vector<Scalar>& upper,
                                   int geqInd,
                                   int leqInd,
                                   std::vector<RangeTreeNode<Scalar,Value>* >& nodes,
                                   std::vector<std::pair<int,int> >& inds) {
            if (leqInd < geqInd) {
                return;
            }

            int compareInd = point->dim() - 2;

            if ((*point)[compareInd] <= upper[compareInd]) {
                if (isLeaf) {
                    nodes.emplace_back(this);
                    inds.emplace_back(std::pair<int,int>(0,0));
                    return;
                }

                int geqIndLeft = pointerToGeqLeft[geqInd];
                int leqIndLeft = pointerToLeqLeft[leqInd];
                if (leqIndLeft >= geqIndLeft) {
                    nodes.emplace_back(left.get());
                    if (left->isLeaf) {
                        inds.emplace_back(std::pair<int,int>(0,0));
                    }
                    else {
                        inds.emplace_back(std::pair<int,int>(geqIndLeft,leqIndLeft));
                    }
                }
                right->rightFractionalCascade(upper,
                                              pointerToGeqRight[geqInd],
                                              pointerToLeqRight[leqInd],
                                              nodes,
                                              inds);
            }
            else {
                if (isLeaf) {
                    return;
                }
                left->rightFractionalCascade(upper,
                                             pointerToGeqLeft[geqInd],
                                             pointerToLeqLeft[leqInd],
                                             nodes,
                                             inds);
            }
        }

        /**
        * Helper function for countInRange(...).
        * @param lower
        * @param withLower
        * @param nodes
        */
        void leftCanonicalNodes(const std::vector<Scalar>& lower,
            std::vector<std::shared_ptr<RangeTreeNode<Scalar, Value> > >& nodes)
        {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (lower[compareInd] <= (*point)[compareInd]) {
                nodes.emplace_back(right);
                if (left->isLeaf) {
                    nodes.emplace_back(left);
                }
                else {
                    left->leftCanonicalNodes(lower, nodes);
                }
            }
            else {
                if (right->isLeaf) {
                    nodes.emplace_back(right);
                }
                else {
                    right->leftCanonicalNodes(lower, nodes);
                }
            }
        }

        /**
        * Helper function for countInRange(...).
        * @param upper
        * @param nodes
        */
        void rightCanonicalNodes(const std::vector<Scalar>& upper,
            std::vector<std::shared_ptr<RangeTreeNode<Scalar, Value> > >& nodes)
        {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (upper[compareInd] >= (*point)[compareInd])
            {
                nodes.emplace_back(left);
                if (right->isLeaf)
                {
                    nodes.emplace_back(right);
                }
                else
                {
                    right->rightCanonicalNodes(upper, nodes);
                }
            }
            else {
                if (left->isLeaf)
                {
                    nodes.emplace_back(left);
                }
                else {
                    left->rightCanonicalNodes(upper, nodes);
                }
            }
        }

        /**
        * Print the structure of the tree rooted at the curret node.
        *
        * The printed structure does not reflect any subtrees for other coordinates.
        *
        * @param numIndents the number of indents to use before every line printed.
        */
        void print(int numIndents) {
            for (int i = 0; i < numIndents; i++) { std::cout << "\t"; }
            if (isLeaf) {
                point->print(true);
            } else {
                point->print(false);
                left->print(numIndents + 1);
                right->print(numIndents + 1);
            }
        }
    };

    /**
    * A class facilitating fast orthogonal range queries.
    *
    * A RangeTree allows for 'orthogonal range queries.' That is, given a collection of
    * points P = {p_1, ..., p_n} in euclidean d-dimensional space, a RangeTree can efficiently
    * answer questions of the form
    *
    * "How many points of p are in the box high dimensional rectangle
    * [l_1, u_1] x [l_2, u_2] x ... x [l_d, u_d]
    * where l_1 <= u_1, ..., l_n <= u_n?"
    *
    * It returns the number of such points in worst case
    * O(log(n)^d) time. It can also return the points that are in the rectangle in worst case
    * O(log(n)^d + k) time where k is the number of points that lie in the rectangle.
    *
    * The particular algorithm implemented here is described in Chapter 5 of the book
    *
    * Mark de Berg, Otfried Cheong, Marc van Kreveld, and Mark Overmars. 2008.
    * Computational Geometry: Algorithms and Applications (3rd ed. ed.). TELOS, Santa Clara, CA, USA.
    */
    template<typename Scalar = double, typename Value = size_t>
    class RangeTree {
        static_assert(std::is_arithmetic<Scalar>::value, "Type Scalar must be numeric");
    private:
        std::shared_ptr<RangeTreeNode<Scalar,Value> > root;
        std::vector<std::shared_ptr<Point<Scalar,Value> > > savedPoints;
        std::vector<Point<Scalar,Value>* > savedPointsRaw;

        std::vector<std::shared_ptr<Point<Scalar,Value> > > copyPointsToHeap(const std::vector<Point<Scalar,Value> >& points) {
            std::vector<std::shared_ptr<Point<Scalar,Value> > > vecOfPointers;
            vecOfPointers.reserve(points.size());
            for (size_t i = 0, length = points.size(); i < length; ++i) {
                vecOfPointers.emplace_back(std::shared_ptr<Point<Scalar,Value> >(new Point<Scalar,Value>(points[i])));
            }
            return vecOfPointers;
        }

        std::vector<Point<Scalar,Value>* > getRawPointers(std::vector<std::shared_ptr<Point<Scalar,Value> > >& points) {
            std::vector<Point<Scalar,Value>* > vecOfPointers;
            vecOfPointers.reserve(points.size());
            for (size_t i = 0, length = points.size(); i < length; ++i) {
                vecOfPointers.emplace_back(points[i].get());
            }
            return vecOfPointers;
        }

        std::vector<Scalar> getModifiedLower(const std::vector<Scalar>& lower,
                         const std::vector<bool>& withLower) const {
            std::vector<Scalar> newLower = lower;
            for (size_t i = 0, length = lower.size(); i < length; ++i) {
                if (std::is_integral<Scalar>::value) {
                    if (!withLower[i]) {
                        ++newLower[i];
                    }
                    }
                else {
                    if (!withLower[i]) {
#if _MSC_VER < 1800
                        newLower[i] = _nextafter(newLower[i], std::numeric_limits<Scalar>::max());
#else
                        newLower[i] = std::nextafter(newLower[i], std::numeric_limits<Scalar>::max());
#endif // 
                    }
                }
            }
            return newLower;
        }

        std::vector<Scalar> getModifiedUpper(const std::vector<Scalar>& upper,
                                        const std::vector<bool>& withUpper) const {
            std::vector<Scalar> newUpper = upper;
            for (size_t i = 0, length = upper.size(); i < length; ++i) {
                if (std::is_integral<Scalar>::value) {
                    if (!withUpper[i]) {
                        --newUpper[i];
                    }
                }
                else {
                    if (!withUpper[i]) {
#if _MSC_VER < 1800
                        newUpper[i] = _nextafter(newUpper[i], std::numeric_limits<Scalar>::lowest());
#else
                        newUpper[i] = std::nextafter(newUpper[i], std::numeric_limits<Scalar>::lowest());
#endif // 
                    }
                }
            }
            return newUpper;
        }

    public:
        /**
        * Construct a new RangeTree from input points.
        *
        * Input points may have duplicates but if two points p_1 and p_2 are at the same euclidean position
        * then they are required to have the same value as all duplicate points will be accumulated into a
        * single point with multiplicity/count equal to the sum of the multiplicities/counts of all such duplicates.
        *
        * @param points the points from which to create a RangeTree
        */
        RangeTree(const std::vector<Point<Scalar,Value> >& points): savedPoints(copyPointsToHeap(points)),
                                                           savedPointsRaw(getRawPointers(savedPoints)) {
            SortedPointMatrix<Scalar,Value> spm(savedPointsRaw);
            root = std::shared_ptr<RangeTreeNode<Scalar,Value> >(new RangeTreeNode<Scalar,Value>(spm));
        }

        /**
        * The number of points within a high dimensional rectangle.
        *
        * The rectangle is defined by the input parameters. In particular, an n-dimensional point
        * p_1 = (p_{11},...,p_{1n}) is an is in the rectangle if, for all 1 <= i <= n, we have
        *
        * lower[i] <= p_{1i} <= upper[i] if withLower[i] == true and withUpper[i] == true, or
        * lower[i] < p_{1i} <= upper[i] if withLower[i] == false and withUpper[i] == true, or
        * lower[i] <= p_{1i} < upper[i] if withLower[i] == true and withUpper[i] == false, or
        * lower[i] < p_{1i} < upper[i] if withLower[i] == false and withUpper[i] == false.
        *
        * @param lower the lower bounds of the rectangle.
        * @param upper the upper bounds of the rectangle.
        * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
        *                  for the lower bounds.
        * @param withUpper as for \withLower but for the upper bounds.
        * @return the number of points in the rectangle.
        */
        int countInRange(const std::vector<Scalar>& lower,
                         const std::vector<Scalar>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            if (lower.size() != upper.size() || lower.size() != withLower.size() ||
                    lower.size() != withUpper.size()) {
                throw std::logic_error("All vectors inputted to countInRange must have the same length.");
            }
            for (size_t i = 0, length = lower.size(); i < length; ++i) {
                if (((!withUpper[i] || !withLower[i]) && lower[i] >= upper[i]) ||
                    lower[i] > upper[i]) {
                    return 0;
                }
            }
            return root->countInRange(getModifiedLower(lower, withLower),
                                      getModifiedUpper(upper, withUpper));
        }

        /**
        * The number of points within a high dimensional rectangle.
        *
        * The rectangle is defined by the input parameters. In particular, an n-dimensional point
        * p_1 = (p_{11},...,p_{1n}) is an is in the rectangle if, for all 1 <= i <= n, we have
        *
        * lower[i] <= p_{1i} <= upper[i]
        *
        * @param lower the lower bounds of the rectangle.
        * @param upper the upper bounds of the rectangle.

        * @return the number of points in the rectangle.
        */
        int countInRange(const std::vector<Scalar>& lower,
                         const std::vector<Scalar>& upper) const {
            if (lower.size() != upper.size()) {
                throw std::logic_error("upper and lower in countInRange must have the same length.");
            }
            return root->countInRange(lower, upper);
        }

        /**
        * Return all points in range.
        *
        * Returns a std::vector of all points in the given rectangle. See \countInRange for how
        * this rectangle is specified. NOTE: these points may not be identical to those points
        * that were given as input to the RangeTree at construction time. This is because
        * duplicate points are merged together with appropriate incrementing of their multiplicity.
        * That is, two points at euclidean position (1,2,3) and multiplicities/counts of 2 and 3
        * respectively will be merged into a single Point with position (1,2,3) and multiplicity 5
        * (recall that all points with the same euclidean position are required to have the same
        * associated value so it is ok to merge in this way).
        *
        * @param lower the lower bounds of the rectangle.
        * @param upper the upper bounds of the rectangle.
        * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
        *                  for the lower bounds.
        * @param withUpper as for \withLower but for the upper bounds.
        * @return the number of points in the rectangle.
        */
        std::vector<Point<Scalar,Value> > pointsInRange(const std::vector<Scalar>& lower,
                                               const std::vector<Scalar>& upper,
                                               const std::vector<bool>& withLower,
                                               const std::vector<bool>& withUpper) const {
            if (lower.size() != upper.size() || lower.size() != withLower.size() ||
                lower.size() != withUpper.size()) {
                throw std::logic_error("All vectors inputted to pointsInRange must have the same length.");
            }
            for (size_t i = 0, length = lower.size(); i < length; ++i) {
                if (((!withUpper[i] || !withLower[i]) && lower[i] >= upper[i]) ||
                    lower[i] > upper[i]) {
                    return std::vector<Point<Scalar,Value> >();
                }
            }
            return root->pointsInRange(getModifiedLower(lower, withLower),
                                       getModifiedUpper(upper, withUpper));
        }

        void print() const {
            root->print(0);
        }
    };

/**
* A class which is used to naively count the number of points in a given rectangle. This class is used
* for testing an benchmarking, it should not be used in practice.
*/
    template<typename Scalar = double, typename Value = size_t>
    class NaiveRangeCounter {
        static_assert(std::is_arithmetic<Scalar>::value, "Type Scalar must be numeric");
    private:
        std::vector<Point<Scalar,Value> > points;

        static bool pointInRange(const Point<Scalar,Value>& point,
                                 const std::vector<Scalar>& lower,
                                 const std::vector<Scalar>& upper,
                                 const std::vector<bool>& withLower,
                                 const std::vector<bool>& withUpper) {
            for (int i = 0; i < point.dim(); ++i) {
                if (point[i] < lower[i] ||
                    (point[i] == lower[i] && !withLower[i])) {
                    return false;
                }
                if (point[i] > upper[i] ||
                    (point[i] == upper[i] && !withUpper[i])) {
                    return false;
                }
            }
            return true;
        }

    public:
        NaiveRangeCounter(std::vector<Point<Scalar,Value> > points): points(points) {}

        int countInRange(const std::vector<Scalar>& lower,
                         const std::vector<Scalar>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            int count = 0;
            for (size_t i = 0, length = points.size(); i < length; ++i) {
                if (pointInRange(points[i], lower, upper, withLower, withUpper)) {
                    count += points[i].count();
                }
            }
            return count;
        }

        std::vector<Point<Scalar,Value> > pointsInRange(const std::vector<Scalar>& lower,
                                               const std::vector<Scalar>& upper,
                                               const std::vector<bool>& withLower,
                                               const std::vector<bool>& withUpper) const {
            std::vector<Point<Scalar,Value> > selectedPoints = {};
            for (size_t i = 0, length = points.size(); i < length; ++i) {
                if (pointInRange(points[i], lower, upper, withLower, withUpper)) {
                    selectedPoints.emplace_back(points[i]);
                }
            }
            return selectedPoints;
        }
    };

} // namespace

#endif //RANGETREE_H
