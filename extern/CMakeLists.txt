include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG ca5835b9d7f9fc544a342eff475a9ba645410ab1
)

FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  manticore
  GIT_REPOSITORY https://github.com/tatami-inc/manticore
  GIT_TAG 2b7e39b39583a2d6bc52b2b4548a06698c9198a2
)

FetchContent_MakeAvailable(manticore)
