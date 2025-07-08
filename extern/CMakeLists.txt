include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master # ^4.0.0
)

FetchContent_Declare(
  tatami_chunked
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_chunked
  GIT_TAG master # ^2.0.0
)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG master # ^1.0.0
)

FetchContent_Declare(
  sanisizer 
  GIT_REPOSITORY https://github.com/LTLA/sanisizer
  GIT_TAG master # ^1.0.0
)

FetchContent_MakeAvailable(manticore)
FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(tatami_chunked)
FetchContent_MakeAvailable(sanisizer)
