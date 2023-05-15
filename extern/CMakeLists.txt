include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/LTLA/tatami
  GIT_TAG ca5835b9d7f9fc544a342eff475a9ba645410ab1
)

FetchContent_MakeAvailable(tatami)
