add_library(Utility
    STATIC
        Utility.c
    )

target_include_directories(Utility
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )