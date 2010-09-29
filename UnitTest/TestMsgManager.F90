program TestMsgManager

    use MsgManager

    call MsgManager_RecordSpeaker("TestMsgManager")
    call MsgManager_Speak(Notice, "Hello, world!")

    call TestSub1
    call TestSub2
    call TestSub3

    call MsgManager_Speak(Notice, "TestSub1 has been removed from speaker stack.")

    call MsgManager_AddConfig("Mod1", "A", "1")
    call MsgManager_AddConfig("Mod1", "B", "2")
    call MsgManager_AddConfig("Mod2", "C", "3")

    call MsgManager_ShowConfig

contains

    subroutine TestSub1

        call MsgManager_RecordSpeaker("TestSub1")
        call MsgManager_Speak(Warning, "I am inside TestMsgManager!")
    
    end subroutine TestSub1

    subroutine TestSub2

        call MsgManager_RecordSpeaker("TestSub2")
        call MsgManager_Speak(Error, "TestSub1 forgot to call delete speaker operation!")
        call MsgManager_DeleteSpeaker
    
    end subroutine TestSub2

    subroutine TestSub3

        call MsgManager_RecordSpeaker("TestSub3")
        call MsgManager_Speak(Notice, "I have fixed it! See below.")
        call MsgManager_DeleteSpeaker
        call MsgManager_DeleteSpeaker
    
    end subroutine TestSub3

end program TestMsgManager
