package be.cmpg.cancer

import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.graph.NetworkReader

object NetworkFactory {
  def loadNetwork(networks: Seq[String], networkFolder: String,excludedGenes:List[String]=List()): Network = {

    val interactionsMap = scala.collection.mutable.Map[Interaction, Set[String]]()

    networks.foreach { n =>
      {
        if (n == "HT") {

          val pendingInteractions = NetworkReader.fromSif(networkFolder + "HumanBinaryHQ_HT.txt", 2, -1, 3, 4)
          pendingInteractions.foreach(interaction => {
            if (interactionsMap.get(interaction).isDefined) {
              interactionsMap.put(interaction, interaction.evidence.++(interactionsMap(interaction)))
            } else (interactionsMap.put(interaction, interaction.evidence))
          })
        } else if (n == "hiII14") {

          val pendingInteractions = NetworkReader.fromTSV(networkFolder + "HI-II-14_Annotated.tsv")._1
          pendingInteractions.foreach(interaction => {
            if (interactionsMap.get(interaction).isDefined) {
              interactionsMap.put(interaction, interaction.evidence.++(interactionsMap(interaction)))
            } else (interactionsMap.put(interaction, interaction.evidence))
          })
        } else if (n == "reactome") {
          // Note this is a very limited and old version of reactome
          val pendingInteractions = NetworkReader.fromSif(networkFolder + "hrn1_reactome.sif", 0, -1, 2, 3)
          pendingInteractions.foreach(interaction => {
            if (interactionsMap.get(interaction).isDefined) {
              interactionsMap.put(interaction, interaction.evidence.++(interactionsMap(interaction)))
            } else (interactionsMap.put(interaction, interaction.evidence))
          })
        } else if (n == "MEMo") {

          val hrn1ppi_i = NetworkReader.fromSif("networks/MEMo/hrn1_ppi.txt", 0, -1, 3)
          val hrn1cellmap_i = NetworkReader.fromSif("networks/MEMo/hrn1_cell-map.sif", 0, -1, 2)
          val hrn1ncinature_i = NetworkReader.fromSif("networks/MEMo/hrn1_nci-nature.sif", 0, -1, 2)
          val hrn1reactome_i = NetworkReader.fromSif("networks/MEMo/hrn1_reactome.sif", 0, -1, 2)
          val pendingInteractions = hrn1ppi_i ++ hrn1cellmap_i ++ hrn1ncinature_i ++ hrn1reactome_i
          pendingInteractions.foreach(interaction => {
            if (interactionsMap.get(interaction).isDefined) {
              interactionsMap.put(interaction, interaction.evidence.++(interactionsMap(interaction)))
            } else (interactionsMap.put(interaction, interaction.evidence))
          })
        }
      }
    }

    val interactions = interactionsMap.map(interaction => {
      val coreInteraction = interaction._1
      //Interaction(coreInteraction.from, coreInteraction.to, coreInteraction.typ, coreInteraction.direction, coreInteraction.regulatory, coreInteraction.probability, evidence = interaction._2.toSet)
      Interaction(coreInteraction.from, coreInteraction.to, coreInteraction.typ, coreInteraction.direction, coreInteraction.regulatory, evidence = interaction._2.toSet)
    }).toSet

    new Network(interactions.filter(interaction =>{!(excludedGenes.contains(interaction.from.name) || excludedGenes.contains(interaction.to.name))}))
  }
}