add_executable(RangeTreeTest)

# 设置公共包含目录
include_directories(${PROJECT_SOURCE_DIR}/../../Doctest/doctest/)

target_include_directories(RangeTreeTest PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    )

target_sources(RangeTreeTest PRIVATE
    Point_test.cpp
    RangeTree_test.cpp
    SortedPointsMatrix_test.cpp
    )