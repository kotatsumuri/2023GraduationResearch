add_library(QD
    STATIC
        qd.cpp
        qd_operator.cpp
    )

target_include_directories(QD
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(util
    STATIC
        util.cpp
    )

target_include_directories(util
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )