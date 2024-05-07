include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master
)

FetchContent_Declare(
  tatami_chunked
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_chunked
  GIT_TAG master
)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG master
)

FetchContent_MakeAvailable(manticore)
FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(tatami_chunked)
