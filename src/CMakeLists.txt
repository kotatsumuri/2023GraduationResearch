add_library(Utility
    STATIC
        Utility.c
    )

target_include_directories(Utility
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(DFT
    STATIC
        DFT.c
    )

target_include_directories(DFT
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )

add_library(FFT
    STATIC
        FFT.c
    )

target_include_directories(FFT
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    )