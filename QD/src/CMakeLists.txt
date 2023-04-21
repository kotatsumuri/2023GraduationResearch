add_library(qd
    STATIC
        qd.cpp
    )

target_include_directories(qd
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(util
    STATIC
        util.cpp
    )

target_include_directories(util
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )