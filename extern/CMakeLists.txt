include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master # ^v3.0.0
)

FetchContent_Declare(
  tatami_chunked
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_chunked
  GIT_TAG master # ^v2.0.0
)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG master # ^v1.0.0
)

FetchContent_Declare(
  subpar 
  GIT_REPOSITORY https://github.com/LTLA/subpar
  GIT_TAG master # ^v0.3.2
)

FetchContent_MakeAvailable(manticore)
FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(tatami_chunked)
FetchContent_MakeAvailable(subpar)
