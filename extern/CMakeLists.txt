include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG 05ba9b011279f900a019cbc255176646bd167b9f
)

FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG 3cf0c49e090179fcb2bbef9b55e011b3bc3c7c85
)

FetchContent_MakeAvailable(manticore)
