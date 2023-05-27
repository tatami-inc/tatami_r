include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG fb428b53b8bfd642db09d71ac34c41351e492b70
)

FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG 2b7e39b39583a2d6bc52b2b4548a06698c9198a2
)

FetchContent_MakeAvailable(manticore)
