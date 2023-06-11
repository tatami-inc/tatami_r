include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG b9f983fcc3c78274aa1404ae554c0c2ae5dedcfc
)

FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG 3cf0c49e090179fcb2bbef9b55e011b3bc3c7c85
)

FetchContent_MakeAvailable(manticore)
