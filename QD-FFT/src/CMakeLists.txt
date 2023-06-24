add_library(util_calc
    STATIC
        util_calc.cpp
    )

target_include_directories(util_calc
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )