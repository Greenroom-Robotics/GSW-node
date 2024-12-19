//
// Created by blake on 9/26/24.
//
#include <BoundTeosBase.h>
#include <BoundTeosIce.h>
#include <BoundTeosSea.h>
#include <napi.h>

Napi::Object Init(Napi::Env env, Napi::Object exports)
{
    Napi::HandleScope scope(env);
    BoundTeosBase::Init(env, exports);
    BoundTeosIce::Init(env, exports);
    BoundTeosSea::Init(env, exports);

    return exports;
}

NODE_API_MODULE(NODE_GYP_MODULE_NAME, Init)
