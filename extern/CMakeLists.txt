include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master
)

FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG master
)

FetchContent_MakeAvailable(manticore)
