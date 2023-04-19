add_library(qd
    STATIC
        qd.cpp
    )

target_include_directories(qd
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(utility
    STATIC
        utility.cpp
    )

target_include_directories(utility
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )