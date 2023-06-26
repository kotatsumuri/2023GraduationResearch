add_library(util_calc
    STATIC
        util_calc.cpp
    )

target_include_directories(util_calc
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(util_io
    STATIC
        util_io.cpp
    )

target_include_directories(util_io
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(qd
    STATIC
        qd.cpp
        qd_util.cpp
    )

target_include_directories(qd
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(fft
    STATIC
        fft.cpp
    )

target_include_directories(fft
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )