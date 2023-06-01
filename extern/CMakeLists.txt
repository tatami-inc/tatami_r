include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG 898039be5c0f9a4ff5de6f95ba7bdfd8d3940119
)

FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG 313002dcea65ddb0fd334034f0f618626c736037
)

FetchContent_MakeAvailable(manticore)
