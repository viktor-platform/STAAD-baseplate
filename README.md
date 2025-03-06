Sub AddRigidControlDependentRelation()
    Dim objOpenSTAAD As Object
    Set objOpenSTAAD = GetObject(, "StaadPro.OpenSTAAD")
    
    Dim retVal As Long
    Dim controlNode As Long
    Dim dependentNodes() As Long ' Use Long array instead of Variant
    
    ' Set the control node (update this number if necessary)
    controlNode = 1414
    
    ' Define and populate the Long array of dependent nodes
    ReDim dependentNodes(0 To 3)
    dependentNodes(0) = 21
    dependentNodes(1) = 399
    dependentNodes(2) = 420
    dependentNodes(3) = 441
    
    ' Call AddControlDependentRelation with corrected syntax
    retVal = objOpenSTAAD.Property.AddControlDependentRelation(controlNode, -1, 1, 1, 1, 1, 1, 1, dependentNodes)
    
    ' Check the return value
    If retVal = 0 Then
        MsgBox "Control/Dependent relation added successfully."
    Else
        MsgBox "Error adding control/dependent relation: " & retVal
    End If
End Sub